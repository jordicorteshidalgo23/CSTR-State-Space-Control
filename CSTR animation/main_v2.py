from manim import *
import numpy as np
import re
from pathlib import Path

# ==========================================
# CONFIGURACIÓN
# ==========================================
DATA_FILE = "simulink_data.csv"   # ver load_simulation() para el formato
SPIN_SPEED = 5.0                  # rad/s de la hélice (ponlo a 0 si quieres el reactor totalmente quieto)
FLOW_SPEED = 0.3                  # vueltas/s de las bolitas de entrada/salida (a caudal máximo)

# Secuencia de la animación (segundos)
DRAW_TIME = 1.2                   # dibujo rápido del reactor
FILL_TIME = 0.8                   # aparece líquido, camisa y anotaciones
HOLD_TIME = 1.0                  # el reactor se queda quieto
FADE_TIME = 0.8                   # desvanecimiento del reactor
AFTER_PLOTS = 0.5                 # pausa con las gráficas dibujadas, antes de que corra la simulación
END_HOLD = 3.0                    # las gráficas se mantienen al terminar la simulación
BLANK_TIME = 0.4                  # pantalla vacía entre el reactor y las gráficas
PLOTS_IN_TIME = 2.5               # las gráficas se dibujan trazo a trazo

VESSEL_W, VESSEL_H = 2.2, 3.2
CORNER_R = 0.3
SURFACE_OFFSET = 0.4

SIG_COLOR = {"cA": BLUE, "cB": RED, "Tc": "#EDF110", "QinV": GREEN_C, "QH": TEAL_C,
             "cAin": "#2BF110", "Tin": PINK, "cBhat": "#2BF110"}

FORMAT = "landscape"

LAYOUTS = {
    "landscape": dict(
        pixels=(1920, 1080), frame=(14.2222, 8.0), fps=60, run_time=12,
        center=(0.0, 0.0), scale=1.3,               # reactor centrado y más grande (ya no comparte pantalla)
        lf=1.0, lf_small=1.0,
        mix_two_lines=False, mix_dx=2.9, out_below=False,
        # Gráficas: 3 paneles apilados (cB, Tc, Tin)
        frame_l=-5.0, frame_w=10.6, frame_h=1.45, tops=(2.8, 0.45, -1.9),
        leg_fs=18, tick_fs=16, ylab_fs=18, title_fs=24, xnum_dy=0.17, yticks=3,
        lw=3.5, frame_lw=1.5,
        xunit_pos=(6.2, -3.52),
    ),
}
L = LAYOUTS[FORMAT]
config.pixel_width, config.pixel_height = L["pixels"]
config.frame_width, config.frame_height = L["frame"]
if L["fps"]:
    config.frame_rate = L["fps"]
RUN_TIME = L["run_time"]


# ==========================================
# MODELO DEL PAPER (Klatt-Engell CSTR) — solo para los datos de DEMO
# ==========================================
V, CP, RHO = 0.01, 3.01, 934.2
DH_AB, DH_BC, DH_AD = 4.2, -11.0, -41.85
AR, KW, MK, CPJ = 0.215, 4032.0, 5.0, 2.0
K10, K20 = 1.287e12, 9.043e6
E1, E2 = -9758.3, -8560.0


def rate_constants(T):
    TK = T + 273.15
    return K10 * np.exp(E1 / TK), K20 * np.exp(E2 / TK)


def plant(x, qv, qh, cAin, Tin):
    cA, cB, T, Tc = x
    k1, k2 = rate_constants(T)
    dcA = qv * (cAin - cA) - k1 * cA - k2 * cA**2
    dcB = -qv * cB + k1 * (cA - cB)
    dT = (qv * (Tin - T)
          - (k1 * cA * DH_AB + k1 * cB * DH_BC + k2 * cA**2 * DH_AD) / (RHO * CP)
          + (Tc - T) * AR * KW / (RHO * CP * V))
    dTc = (qh + AR * KW * (T - Tc)) / (MK * CPJ)
    return np.array([dcA, dcB, dT, dTc])


def smooth_step(t, t_step, width=0.02):
    return 0.5 * (1 + np.tanh((t - t_step) / width))


def demo_inputs(t):
    qv = NOMINAL["QinV"] + (22.0 - NOMINAL["QinV"]) * smooth_step(t, 0.4)
    qh = NOMINAL["QH"] + (-1800.0 - NOMINAL["QH"]) * smooth_step(t, 0.8)
    cAin = NOMINAL["cAin"] + 100.0 * np.sin(8 * t) * smooth_step(t, 1.2)
    Tin = NOMINAL["Tin"] + 3.0 * smooth_step(t, 1.5)
    return qv, qh, cAin, Tin


def demo_simulation(t_end=2.0, dt=0.0005):
    x = np.array([NOMINAL[k] for k in ("cA", "cB", "T", "Tc")])
    ts = np.linspace(0, t_end, int(round(t_end / dt)) + 1)
    X = np.empty((len(ts), 4))
    X[0] = x
    for i in range(len(ts) - 1):
        t = ts[i]
        f = lambda tt, xx: plant(xx, *demo_inputs(tt))
        k1 = f(t, x); k2 = f(t + dt / 2, x + dt / 2 * k1)
        k3 = f(t + dt / 2, x + dt / 2 * k2); k4 = f(t + dt, x + dt * k3)
        x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        X[i + 1] = x
    qv, qh, cAin, Tin = (np.array(v, dtype=float) * np.ones_like(ts)
                         for v in zip(*[demo_inputs(t) for t in ts]))
    err = np.exp(-ts / 0.05)
    return SimData(ts, cA=X[:, 0], cB=X[:, 1], T=X[:, 2], Tc=X[:, 3],
                   QinV=qv, QH=qh, cAin=cAin, Tin=Tin,
                   cAhat=X[:, 0] - 400.0 * err, cBhat=X[:, 1] + 250.0 * err)


# ==========================================
# DATOS (SIMULINK O DEMO)
# ==========================================
SPECIES = ["cA", "cB", "cC", "cD"]
SIGNALS = SPECIES + ["T", "Tc", "QinV", "QH", "cAin", "Tin", "cAhat", "cBhat"]
NOMINAL = {"cA": 2140.0, "cB": 1090.0, "T": 114.2, "Tc": 112.9,
           "QinV": 14.19, "QH": -1113.5, "cAin": 5100.0, "Tin": 104.9}


def _byproducts(t, cA, cB, T, q):
    k1, k2 = rate_constants(T)
    src = [k1 * cB, 0.5 * k2 * cA**2]
    out = []
    for s in src:
        c = np.zeros_like(t)
        c[0] = s[0] / max(q[0], 1e-9)
        for i in range(len(t) - 1):
            qi = max(q[i], 1e-9)
            e = np.exp(-qi * (t[i + 1] - t[i]))
            c[i + 1] = c[i] * e + s[i] / qi * (1 - e)
        out.append(c)
    return out


class SimData:
    """Señales vs tiempo. Acceso: sim.at('QinV', t), sim.conc(t) -> (...,4) [A,B,C,D]."""

    def __init__(self, t, **signals):
        t = np.asarray(t, float)
        order = np.argsort(t, kind="stable")
        self.t = t[order]
        self.s = {}
        for k, v in signals.items():
            self.s[k] = np.asarray(v, float)[order]
        for k in ("cA", "cB"):
            if k not in self.s:
                raise ValueError(f"Falta la columna obligatoria '{k}'")
        for k, v in NOMINAL.items():
            self.s.setdefault(k, np.full_like(self.t, v))
        if "cC" not in self.s or "cD" not in self.s:
            self.s["cC"], self.s["cD"] = _byproducts(
                self.t, self.s["cA"], self.s["cB"], self.s["T"], self.s["QinV"])
        for k in SPECIES:
            self.s[k] = np.maximum(self.s[k], 0)
        self.t0, self.tf = self.t[0], self.t[-1]
        self.cmax = float(max(self.s[k].max() for k in SPECIES))

    def at(self, name, t):
        return np.interp(t, self.t, self.s[name])

    def conc(self, t):
        return np.stack([self.at(k, t) for k in SPECIES], axis=-1)


_EXACT = {"t"} | set(SIGNALS)
_LOOSE = {re.sub(r"[^a-z0-9]", "", k.lower()): k for k in SIGNALS if k != "T"}
_LOOSE.update(time="t", tout="t", caest="cAhat", cbest="cBhat")


def _canon(name):
    if name in _EXACT:
        return name
    return _LOOSE.get(re.sub(r"[^a-z0-9]", "", name.lower()))


def load_simulation(path):
    """
    CSV con cabecera (el orden no importa). Obligatorias: t, cA, cB.
    Opcionales: T, Tc, QinV, QH, cAin, Tin, cAhat, cBhat, cC, cD.
    Unidades: t [hr], c [mol/m3], T [°C], Qin/V [1/hr], QH [kJ/hr].
    """
    p = Path(path)
    raw = {}
    if p.suffix == ".mat":
        from scipy.io import loadmat
        items = {k: v for k, v in loadmat(p).items() if not k.startswith("__")}
    else:
        data = np.genfromtxt(p, delimiter=",", names=True, encoding=None)
        items = {n: data[n] for n in data.dtype.names}
    for name, val in items.items():
        k = _canon(name)
        if k:
            raw[k] = np.asarray(val, float).squeeze()
    if "t" not in raw:
        raise ValueError("Falta la columna de tiempo 't'")
    return SimData(raw.pop("t"), **raw)


# ==========================================
# UTILIDADES DE EJES
# ==========================================
def nice_step(span, n=5):
    raw = span / n
    mag = 10 ** np.floor(np.log10(raw))
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return float(m * mag)


def ticks(lo, hi, step):
    dec = max(6, int(-np.floor(np.log10(step))) + 3)
    n = int(np.floor((hi - lo) / step + 1e-9))
    vals = np.round(lo + step * np.arange(n + 1), dec)
    places = int(max(0, -np.floor(np.log10(step) + 1e-9)))
    return [float(v) for v in vals], places


def nice_range(arrays, n=4, pad=0.08, min_free=0.0):
    lo = min(float(np.min(a)) for a in arrays)
    hi = max(float(np.max(a)) for a in arrays)
    dmax = hi
    if hi - lo < 1e-9 * max(1.0, abs(hi)):
        d = max(abs(hi) * 0.05, 0.5)
        lo, hi = lo - d, hi + d
    else:
        d = (hi - lo) * pad
        lo, hi = lo - d, hi + d
    step = nice_step(hi - lo, n)
    dec = max(6, int(-np.floor(np.log10(step))) + 3)
    lo_n = round(float(np.floor(lo / step) * step), dec)
    hi_n = round(float(np.ceil(hi / step) * step), dec)
    for _ in range(20):
        if (hi_n - dmax) >= min_free * (hi_n - lo_n):
            break
        hi_n = round(hi_n + step, dec)
    return lo_n, hi_n, step


class JacketedCSTR(Scene):
    def construct(self):
        # ---------- Datos ----------
        sim = load_simulation(DATA_FILE) if Path(DATA_FILE).exists() else demo_simulation()
        tracker = ValueTracker(sim.t0)
        now = tracker.get_value

        qv_ref = max(float(np.abs(sim.s["QinV"]).max()), 1e-9)
        qh_ref = max(float(np.abs(sim.s["QH"]).max()), 1e-9)
        temps = np.concatenate([sim.s["T"], sim.s["Tc"], sim.s["Tin"]])
        T_lo, T_hi = float(temps.min()), float(temps.max())
        if T_hi - T_lo < 1e-6:
            T_hi = T_lo + 1.0

        def temp_color(T):
            f = float(np.clip((T - T_lo) / (T_hi - T_lo), 0, 1))
            return interpolate_color(BLUE_C, WHITE, 2 * f) if f < 0.5 else interpolate_color(WHITE, RED_C, 2 * f - 1)

        def heat_color(qh):
            s = float(np.clip(abs(qh) / qh_ref, 0, 1))
            return interpolate_color(WHITE, TEAL_C if qh < 0 else ORANGE, s)

        # ---------- Geometría ----------
        S = L["scale"]
        center = np.array([L["center"][0], L["center"][1], 0.0])
        o = lambda x, y: center + S * np.array([x, y, 0.0])
        fs = lambda base, k="lf": int(round(base * L[k]))
        corner_r = CORNER_R * S

        inner_vessel = RoundedRectangle(corner_radius=corner_r, height=VESSEL_H * S, width=VESSEL_W * S,
                                        color=WHITE, stroke_width=3.5).move_to(center)
        outer_jacket = RoundedRectangle(corner_radius=0.45 * S, height=3.7 * S, width=3.2 * S,
                                        color=WHITE, stroke_width=3.5).move_to(center)

        inlet_path = VMobject(color=WHITE, stroke_width=3.5).set_points_as_corners(
            [o(-1.7, 2.3), o(-0.6, 2.3), o(-0.6, 1.6)])
        inlet_flange = Line(o(-1.7, 2.45), o(-1.7, 2.15), color=WHITE, stroke_width=4)
        outlet_path = VMobject(color=WHITE, stroke_width=3.5).set_points_as_corners(
            [o(0.6, -1.6), o(0.6, -2.3), o(1.7, -2.3)])
        outlet_flange = Line(o(1.7, -2.45), o(1.7, -2.15), color=WHITE, stroke_width=4)

        x_out, x_in = o(-2.15, 0)[0], o(-1.25, 0)[0]
        y_start, y_end = o(0, 0.6)[1], o(0, -0.6)[1]
        n_arcs = 4
        dy = (y_start - y_end) / n_arcs
        top_pipe = Line([x_out, y_start, 0], [x_in, y_start, 0], color=WHITE, stroke_width=3)
        bot_pipe = Line([x_in, y_end, 0], [x_out, y_end, 0], color=WHITE, stroke_width=3)
        arcs = VGroup(*[
            ArcBetweenPoints([x_in, y_start - i * dy, 0], [x_in, y_start - (i + 1) * dy, 0],
                             angle=PI, color=WHITE, stroke_width=3)
            for i in range(n_arcs)
        ])
        coil = VGroup(top_pipe, arcs, bot_pipe)

        gland = Rectangle(height=0.25 * S, width=0.15 * S, color=WHITE, fill_color=BLACK,
                          fill_opacity=1, stroke_width=3).move_to(o(0, 1.85))
        shaft = Line(o(0, 2.5), o(0, -0.8), color=WHITE, stroke_width=4)
        hub = shaft.get_end()

        blade_left = Ellipse(width=0.55 * S, height=0.18 * S, color=WHITE, stroke_width=3).move_to(hub + LEFT * 0.28 * S)
        blade_right = Ellipse(width=0.55 * S, height=0.18 * S, color=WHITE, stroke_width=3).move_to(hub + RIGHT * 0.28 * S)
        propeller = VGroup(blade_left, blade_right)

        def spin(mob, dt):
            mob.rotate(dt * SPIN_SPEED, axis=UP, about_point=hub)

        propeller.add_updater(spin)

        top_spin_arrow = CurvedArrow(o(-0.25, 2.35), o(0.25, 2.35),
                                     radius=0.3 * S, color=WHITE, tip_length=0.15 * S)

        reactor_rig = VGroup(outer_jacket, inner_vessel, inlet_path, inlet_flange,
                             outlet_path, outlet_flange, coil, gland, shaft,
                             propeller, top_spin_arrow)

        # ---------- Efectos (el reactor se muestra en el estado inicial t0) ----------
        jacket_fill = Cutout(outer_jacket.copy(), inner_vessel.copy(),
                             fill_color=WHITE, fill_opacity=0.2, stroke_width=0)
        jacket_fill.add_updater(lambda m: m.set_fill(temp_color(sim.at("Tc", now())), opacity=0.2))
        outer_jacket.add_updater(lambda m: m.set_stroke(temp_color(sim.at("Tc", now()))))
        inner_vessel.add_updater(lambda m: m.set_stroke(temp_color(sim.at("T", now()))))

        def coil_upd(m):
            qh = float(sim.at("QH", now()))
            m.set_stroke(heat_color(qh), width=3 + 3 * min(abs(qh) / qh_ref, 1))

        coil.add_updater(coil_upd)

        def make_flow(path, temp_name, n=5):
            dots = VGroup(*[Dot(radius=0.045) for _ in range(n)])
            dots.phase = 0.0

            def upd(m, dt):
                speed = FLOW_SPEED * float(sim.at("QinV", now())) / qv_ref
                m.phase = (m.phase + dt * speed) % 1
                col = temp_color(sim.at(temp_name, now()))
                for i, d in enumerate(m):
                    d.move_to(path.point_from_proportion((m.phase + i / n) % 1))
                    d.set_color(col)

            upd(dots, 0)
            dots.add_updater(upd)
            return dots

        inlet_dots = make_flow(inlet_path, "Tin")
        outlet_dots = make_flow(outlet_path, "T")

        # ---------- Anotaciones ----------
        in_arrow = Arrow(o(-2.5, 2.3), o(-1.8, 2.3), color=WHITE, buff=0)
        lbl_in_A = MathTex("A", font_size=fs(28)).next_to(in_arrow, UP, buff=0.1)
        lbl_in_params = MathTex(r"Q_{in}", r",\,", r"c_{A,in}", r",\,", r"T_{in}",
                                font_size=fs(24)).next_to(in_arrow, DOWN, buff=0.3).shift(LEFT * 0.3 * S)
        out_arrow = Arrow(o(1.8, -2.3), o(2.6, -2.3), color=WHITE, buff=0)
        lbl_out = MathTex("A, B, C, D", font_size=fs(26))
        if L["out_below"]:
            lbl_out.next_to(out_arrow, DOWN, buff=0.1)
        else:
            lbl_out.next_to(out_arrow, RIGHT, buff=0.15)
        lbl_QH = MathTex("Q_H", font_size=fs(28), color=WHITE).next_to(coil, LEFT, buff=0.2)
        lbl_Tc = MathTex("T_c", font_size=fs(26, "lf_small"), color=WHITE).move_to(o(-1.38, -1.1))
        if L["mix_two_lines"]:
            lbl_mix = VGroup(MathTex(r"c_A, c_B, c_C,", font_size=fs(24)),
                             MathTex(r"c_D, T", font_size=fs(24))).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        else:
            lbl_mix = MathTex(r"c_A, c_B, c_C, c_D, T", font_size=fs(24))
        lbl_mix.move_to(o(L["mix_dx"], 0.3))
        pointer = Line(lbl_mix.get_left() + LEFT * 0.1, o(0.4, 0.05), color=WHITE, stroke_width=2)
        annotations = VGroup(in_arrow, lbl_in_A, lbl_in_params, out_arrow, lbl_out,
                             lbl_QH, lbl_Tc, lbl_mix, pointer)

        # ---------- Líquido ----------
        W = inner_vessel.width
        y_bottom = inner_vessel.get_bottom()[1]
        y_surface = center[1] + SURFACE_OFFSET * S
        H = y_surface - y_bottom
        x_min = inner_vessel.get_left()[0]

        def style(color, fill):
            return VMobject(color=color, fill_color=fill, fill_opacity=0.75, stroke_width=0.5)

        rect_A, rect_B = style(BLUE, BLUE_C), style(GOLD, GOLD_C)
        rect_C, rect_D = style(RED, RED_C), style(PURPLE, PURPLE_C)
        labels = {k: MathTex(k, font_size=fs(22), color=WHITE) for k in "ABCD"}

        x_max = x_min + W

        def y_low(x):
            dl = (x_min + corner_r) - x
            dr = x - (x_max - corner_r)
            dx = max(dl, dr, 0.0)
            if dx <= 0:
                return y_bottom
            return y_bottom + corner_r - np.sqrt(max(corner_r**2 - dx**2, 0.0))

        def clipped_points(x0, x1, y0, y1):
            xs = np.linspace(x0, x1, 60)
            lower = [(x, min(max(y0, y_low(x)), y1)) for x in xs]
            upper = [(x, y1) for x in xs[::-1]]
            pts = [[x, y, 0] for x, y in lower + upper]
            pts.append(pts[0])
            return np.array(pts)

        def update_liquid(mob):
            ca, cb, cc, cd = sim.conc(now())
            total = max(ca + cb + cc + cd, 1e-9)
            h_top, h_bot = H * (ca + cb) / total, H * (cc + cd) / total
            f_ab = ca / (ca + cb) if ca + cb > 1e-6 else 0.5
            f_cd = cc / (cc + cd) if cc + cd > 1e-6 else 0.5
            w_A, w_B = W * f_ab, W * (1 - f_ab)
            w_C, w_D = W * f_cd, W * (1 - f_cd)

            def put(rect, lbl, x0, w, y0, h):
                a, b = x_min + x0, x_min + x0 + w
                c, d = y_bottom + y0, y_bottom + y0 + h
                rect.set_points_as_corners(clipped_points(a, max(b, a + 1e-3), c, max(d, c + 1e-3)))
                lbl.move_to([(a + b) / 2, (c + d) / 2, 0])
                lbl.set_opacity(1 if (w > 0.15 and h > 0.15) else 0)

            put(rect_C, labels["C"], 0, w_C, 0, h_bot)
            put(rect_D, labels["D"], w_C, w_D, 0, h_bot)
            put(rect_A, labels["A"], 0, w_A, h_bot, h_top)
            put(rect_B, labels["B"], w_A, w_B, h_bot, h_top)

        liquid = VGroup(rect_A, rect_B, rect_C, rect_D, *labels.values())
        liquid.add_updater(update_liquid)
        update_liquid(liquid)

        liquid_surface = Line([x_min, y_surface, 0], [x_min + W, y_surface, 0],
                              color=WHITE, stroke_width=2)

        # ---------- Gráficas: 2 paneles ----------
        #   1) "Variables a controlar": eje izq. cB (real, estimada, referencia); eje der. Tc (real, referencia)
        #   2) "Perturbación en Tin": Tin y su referencia (104.9 °C)
        t0, tf = sim.t0, sim.tf
        FRAME_L, FRAME_W = L["frame_l"], L["frame_w"]
        t_plot = np.linspace(t0, tf, 1500)
        x_ticks, x_places = ticks(t0, tf, nice_step(tf - t0, 4))
        px = lambda t: FRAME_L + (np.asarray(t) - t0) / (tf - t0) * FRAME_W
        has_est = "cBhat" in sim.s
        # Referencias constantes = régimen nominal
        sim.s["cBref"] = np.full_like(sim.t, NOMINAL["cB"])
        sim.s["Tcref"] = np.full_like(sim.t, NOMINAL["Tc"])
        sim.s["Tinref"] = np.full_like(sim.t, NOMINAL["Tin"])      # Tin = 104.9 °C

        def make_curve(name, color, style_, py):
            """style_: False (sólida), 'thin' (sólida más fina) o 'ref' (discontinua de referencia)."""
            pts = np.stack([px(t_plot), py(sim.at(name, t_plot)), np.zeros_like(t_plot)], axis=1)
            is_ref = style_ == "ref"
            width = L["lw"] * (0.6 if is_ref else 0.8 if style_ == "thin" else 1.0)
            curve = VMobject(color=color, stroke_width=width, stroke_opacity=0.75 if is_ref else 1.0)
            d_on, d_off = 45, 30                       # solo se usa en las referencias (rectas)

            def upd(m):
                if not is_ref:
                    k = max(int(np.searchsorted(t_plot, now(), side="right")), 2)
                    m.set_points_as_corners(pts[:k])
                    return
                m.clear_points()
                i = 0
                while i < len(pts) - 1:
                    seg = pts[i:min(i + d_on + 1, len(pts))]
                    m.start_new_path(seg[0])
                    for p in seg[1:]:
                        m.add_line_to(p)
                    i += d_on + d_off

            upd(curve)
            curve.add_updater(upd)
            return curve

        def build_panel(top, h, left, right=None, xnums=False, ylabels=(None, None)):
            bottom = top - h
            frame = VGroup(Rectangle(width=FRAME_W, height=h, color=GRAY, stroke_width=L["frame_lw"]
                                     ).move_to([FRAME_L + FRAME_W / 2, top - h / 2, 0]))
            curves = VGroup()
            for xv in x_ticks:
                x = float(px(xv))
                frame.add(Line([x, bottom, 0], [x, top, 0], color=GRAY, stroke_width=1, stroke_opacity=0.15),
                          Line([x, bottom, 0], [x, bottom + 0.06, 0], color=GRAY, stroke_width=1.5))
                if xnums:
                    frame.add(DecimalNumber(xv, num_decimal_places=x_places, font_size=L["tick_fs"],
                                            color=GRAY_B).move_to([x, bottom - L["xnum_dy"], 0]))
            sides = ((left, FRAME_L, LEFT, ylabels[0]), (right, FRAME_L + FRAME_W, RIGHT, ylabels[1]))
            for specs, x0, direction, ylab in sides:
                if not specs:
                    continue
                lo, hi, step = nice_range([sim.at(n, t_plot) for n, _, _ in specs], n=L["yticks"])
                py = lambda v, lo=lo, hi=hi: bottom + (np.asarray(v) - lo) / (hi - lo) * h
                lab_color = specs[0][1] if right else GRAY_B
                vals, places = ticks(lo, hi, step)
                tick_labels = []
                for v in vals:
                    y = float(py(v))
                    lab = DecimalNumber(v, num_decimal_places=places, font_size=L["tick_fs"], color=lab_color
                                        ).next_to([x0, y, 0], direction, buff=0.08)
                    tick_labels.append(lab)
                    frame.add(Line([x0, y, 0], [x0 - direction[0] * 0.06, y, 0], color=GRAY,
                                   stroke_width=1.5), lab)
                if ylab:
                    gap = 0.08 + max(l.width for l in tick_labels) + 0.22
                    frame.add(MathTex(ylab, font_size=L["ylab_fs"], color=lab_color).rotate(PI / 2).move_to(
                        [x0 + direction[0] * gap, top - h / 2, 0]))
                for name, color, st in specs:
                    curves.add(make_curve(name, color, st, py))
            cursor = Line([0, bottom, 0], [0, top, 0], color=WHITE, stroke_width=1, stroke_opacity=0.5)
            cursor.add_updater(lambda m: m.set_x(float(px(now()))))
            return frame, curves, cursor

        def legend_row(items, top, title_txt):
            rows = []
            for color, sym, st in items:
                if st == "ref":
                    sample = DashedLine(ORIGIN, RIGHT * 0.45, color=color, dash_length=0.12,
                                        stroke_width=L["lw"] * 0.6, stroke_opacity=0.75)
                else:
                    sample = Line(ORIGIN, RIGHT * 0.45, color=color,
                                  stroke_width=L["lw"] * (0.8 if st == "thin" else 1.0))
                rows.append(VGroup(sample, MathTex(sym, font_size=L["leg_fs"])).arrange(RIGHT, buff=0.12))
            leg = VGroup(*rows).arrange(RIGHT, buff=0.45)
            cx = FRAME_L + FRAME_W / 2
            leg.move_to([cx, top + 0.1 + leg.height / 2, 0])
            title = Tex(title_txt, font_size=L["title_fs"], color=WHITE)
            title.move_to([cx, leg.get_top()[1] + 0.1 + title.height / 2, 0])
            return VGroup(leg, title)

        TOPS, PH = L["tops"], L["frame_h"]
        C = SIG_COLOR
        items_cB = [(C["cB"], r"c_B\text{ (plant)}", False)] \
            + ([(C["cBhat"], r"c_B\text{ (observer)}", "thin")] if has_est else []) \
            + [(C["cB"], r"c_{B,\text{ref}}", "ref")]
        items_Tc = [(C["Tc"], "T_c", False), (C["Tc"], r"T_{c,\text{ref}}", "ref")]
        items_Tin = [(C["Tin"], r"T_{in}", False), (C["Tin"], r"T_{in,\text{nom}}", "ref")]

        p1 = build_panel(TOPS[0], PH,
                         [("cBref", C["cB"], "ref"), ("cB", C["cB"], False)]
                         + ([("cBhat", C["cBhat"], "thin")] if has_est else []),
                         ylabels=(r"c_B\;(\text{mol/m}^3)", None))
        p2 = build_panel(TOPS[1], PH, [("Tcref", C["Tc"], "ref"), ("Tc", C["Tc"], False)],
                         ylabels=(r"T_c\;(^\circ\text{C})", None))
        p3 = build_panel(TOPS[2], PH, [("Tinref", C["Tin"], "ref"), ("Tin", C["Tin"], False)],
                         xnums=True, ylabels=(r"T_{in}\;(^\circ\text{C})", None))

        heads = VGroup(legend_row(items_cB, TOPS[0], r"Concentration of B ($c_B$)"),
                       legend_row(items_Tc, TOPS[1], r"Jacket temperature ($T_c$)"),
                       legend_row(items_Tin, TOPS[2], r"Disturbance in reactor feed temperature ($T_{in}$)"))
        x_unit = MathTex(r"t\;(\text{hr})", font_size=L["leg_fs"], color=GRAY_B).move_to([*L["xunit_pos"], 0])

        plots = VGroup(*p1, *p2, *p3, heads, x_unit)

        # ---------- Capas ----------
        for mob, z in ((liquid, 1), (liquid_surface, 1), (jacket_fill, 1), (reactor_rig, 2),
                       (inlet_dots, 3), (outlet_dots, 3), (annotations, 3), (plots, 4)):
            mob.set_z_index(z)

        for mob in (outer_jacket, inner_vessel, jacket_fill, coil, liquid, inlet_dots, outlet_dots, plots):
            mob.update(0)

        # ---------- Secuencia ----------
        # 1) el reactor se dibuja rápido
        self.play(Create(reactor_rig, lag_ratio=0.1), run_time=DRAW_TIME)
        # 2) aparecen líquido, camisa, partículas y anotaciones
        self.play(AnimationGroup(FadeIn(jacket_fill), FadeIn(liquid), FadeIn(liquid_surface),
                                 FadeIn(inlet_dots), FadeIn(outlet_dots),
                                 LaggedStart(*[FadeIn(a) for a in annotations], lag_ratio=0.12)),
                  run_time=FILL_TIME)
        # 3) el reactor se queda quieto
        self.wait(HOLD_TIME)
        # 4) el reactor se desvanece entero a la vez y la pantalla queda vacía un instante
        reactor_all = VGroup(reactor_rig, liquid, liquid_surface, jacket_fill,
                             inlet_dots, outlet_dots, annotations)
        reactor_all.clear_updaters()      # si no, los updaters (relleno de camisa, líquido...) pelean con el fundido
        self.play(FadeOut(reactor_all), run_time=FADE_TIME)
        self.clear()
        self.wait(BLANK_TIME)
        # 5) las gráficas se dibujan: marco y ticks, títulos/leyendas y referencias, panel a panel
        panel_anims = []
        for (frame, curves, cursor), head in zip((p1, p2, p3), heads):
            panel_anims.append(AnimationGroup(Create(frame, lag_ratio=0), Write(head),
                                              Create(VGroup(curves, cursor), lag_ratio=0),
                                              lag_ratio=0))
        self.play(LaggedStart(*panel_anims, Write(x_unit), lag_ratio=0.35), run_time=PLOTS_IN_TIME)
        self.wait(AFTER_PLOTS)
        # 6) simulación: las curvas avanzan con el tiempo
        self.play(tracker.animate.set_value(tf), run_time=RUN_TIME, rate_func=linear)
        self.wait(END_HOLD)