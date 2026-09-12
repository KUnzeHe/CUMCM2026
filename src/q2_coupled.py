"""
2026 高教社杯 A题 —— 问题2 求解器（变物性热-湿耦合，固定半径）

模型口径见 PROJECT.md §4.3 与 §5：整个时域统一使用附件3 的 rho(C)、cp(C)、
k(C) 和 D(C, T_K)（扩散系数中的温度用 K），不考虑蒸发潜热；h、hm 沿用附件2。

离散：顶点型守恒有限体积（节点落在 0、0.1、…、2.0 cm，输出免插值）；界面 D 取
调和平均、k 取算术平均；后向欧拉；每个时间层内 Picard 交替更新 C -> 物性 -> T；
三对角方程用 scipy.linalg.solve_banded 直接求解。

入库说明（2026-09-12）：本文件由项目外的 Q2 原型并入，改动仅限于：
  * 附件路径改为 data/official，.npy 输出改为 outputs/q2/data；
  * 函数与参数改名，与 src/q3_drying_time.py 既有的 import 对齐
    （rho_f/cp_f/k_f/D_f/assemble/bsolve、P.T0/P.C0）；
  * 删除潜热分支（eta、Lv、rho_d）。原型中 eta 默认 0，qlat 恒为 0，删除不改变
    任何数值；同时与 §4.3「问题2 不考虑蒸发潜热、不建立含潜热对照」一致；
  * Picard 达到迭代上限时告警，不再静默接受；
  * 输出时刻的匹配改为定点比较，避免 dt != 1 时静默丢输出。
离散格式、时间推进、界面平均方式与物性公式未动。

待办：界面 D 仍为调和平均。按 PROJECT.md §12 第7条，问题2-4 应改用 Kirchhoff
积分平均；该修改同时影响本文件与 q3_drying_time.py，单独作为一步进行。
"""
import os
import warnings

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.linalg import solve_banded

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ENV_PATH = os.path.join(PROJECT_ROOT, 'data', 'official', 'A题',
                                '附件', '附件1.xlsx')
DEFAULT_OUTPUT_DATA = os.path.join(PROJECT_ROOT, 'outputs', 'q2', 'data')


class P:
    R = 0.02            # m          圆柱半径
    L = 0.25            # m          圆柱长度（一维径向模型中不参与计算）
    h = 25.0            # W/(m^2 K)  表面对流换热系数
    hm = 8e-7           # m/s        表面传质系数
    T0 = 28.0           # °C         初始温度
    C0 = 2.55           # kg/kg      初始干基含水率
    t_end = 10800.0     # s          问题2 要求输出前 3h
    N = 200             # 径向控制体数，节点 0..N
    dt = 1.0            # s
    tol = 1e-9          # Picard 收敛判据
    maxit = 40
    # 以下常物性仅供 q1_mode 退化回归检验使用
    Rho = 820.0
    Cp = 2600.0
    k = 0.36


def rho_f(C):
    return 650.0 + 128.0 * C


def cp_f(C):
    return 1450.0 + 2736.0 * C / (C + 1.0)


def k_f(C):
    return 0.21 + 0.38 * C / (C + 1.0)


def D_f(C, T):
    """有效水分扩散系数。T 传入摄氏度，内部换算为 K。"""
    Cs = np.maximum(C, 1e-6)
    Tk = np.maximum(T + 273.15, 1.0)
    return 2.4e-3 * np.exp(-0.45 / Cs) * np.exp(-3850.0 / Tk)


def D_q1(C):
    """问题1 的扩散系数，仅供 q1_mode 退化回归检验使用。"""
    return 7e-9 * np.exp(-0.89 / np.maximum(C, 1e-6))


def load_env(path=None):
    """附件1 分段线性插值；问题2 只用到 0-3h，落在附件范围内，不涉及外推。"""
    if path is None:
        path = DEFAULT_ENV_PATH
    d = pd.read_excel(path).values.astype(float)
    t, T, C = d[:, 0], d[:, 1], d[:, 2]
    return (interp1d(t, T, 'linear', bounds_error=False, fill_value=(T[0], T[-1])),
            interp1d(t, C, 'linear', bounds_error=False, fill_value=(C[0], C[-1])))


def make_grid(N, R):
    """顶点型网格。V 为控制体的面积权重（已除去 2*pi*L）。"""
    dr = R / N
    r = np.arange(N + 1) * dr
    V = r * dr
    V[0] = dr ** 2 / 8                # 中心半控制体
    V[-1] = (dr / 2) * (R - dr / 4)   # 表面半控制体
    rf = (np.arange(N) + 0.5) * dr    # 界面位置
    assert abs(V.sum() - R ** 2 / 2) < 1e-13 * R ** 2, '体积因子之和错误'
    return dr, r, V, rf


def assemble(kf, rf, dr, N, s):
    """守恒型三对角算子。内部行和为 0，表面行和为 -s（s = h*R 或 hm*R）。"""
    A = kf * rf / dr
    di = np.zeros(N + 1)
    di[:-1] -= A
    di[1:] -= A
    di[-1] -= s
    return A.copy(), di, A.copy()


def check_K(lo, di, up, s, N):
    rows = di.copy()
    rows[:-1] += up
    rows[1:] += lo
    assert np.abs(rows[:N]).max() < 1e-9, '内部行和不为 0'
    assert abs(rows[N] + s) < 1e-9, '表面边界系数错误'


def bsolve(lo, di, up, M_dt, rhs, N):
    ab = np.zeros((3, N + 1))
    ab[0, 1:] = -up
    ab[1, :] = M_dt - di
    ab[2, :-1] = -lo
    return solve_banded((1, 1), ab, rhs)


def harm(x):
    return 2 * x[:-1] * x[1:] / (x[:-1] + x[1:])


def arith(x):
    return 0.5 * (x[:-1] + x[1:])


def solve(t_end, T_inf, C_inf, N=None, dt=None, out_times=None,
          q1_mode=False, verbose=False):
    """求解 [0, t_end] 上的耦合温度场与含水率场。

    q1_mode=True 时冻结为问题1 的常物性与常温扩散系数，用于退化回归检验。
    out_times 为 None 时逐步记录；否则只记录集合中的时刻（单位 s）。
    返回 (out, diag)，out[t] = (T, C)。
    """
    N = P.N if N is None else N
    dt = P.dt if dt is None else dt
    dr, r, V, rf = make_grid(N, P.R)
    if out_times is not None:
        out_times = {round(float(x), 6) for x in out_times}

    T = np.full(N + 1, P.T0)
    C = np.full(N + 1, P.C0)
    nstep = int(round(t_end / dt))
    out = {}
    itmax, emax, mmax = 0, 0.0, 0.0
    checked = False
    log_every = max(1, int(round(3600.0 / dt)))

    for n in range(1, nstep + 1):
        t1 = n * dt
        T_old, C_old = T.copy(), C.copy()
        T_new, C_new = T.copy(), C.copy()

        for it in range(P.maxit):
            # (a) 传质：D 在时间中点取值
            if q1_mode:
                Dn = D_q1(0.5 * (C_old + C_new))
            else:
                Dn = D_f(0.5 * (C_old + C_new), 0.5 * (T_old + T_new))
            lo, di, up = assemble(harm(Dn), rf, dr, N, P.hm * P.R)
            if not checked:
                check_K(lo, di, up, P.hm * P.R, N)
            bC = V / dt * C_old
            bC[-1] += P.hm * P.R * C_inf(t1)
            C_next = bsolve(lo, di, up, V / dt, bC, N)

            # (b) 传热：物性用刚更新的 C
            C_mid = 0.5 * (C_old + C_next)
            if q1_mode:
                M = np.full(N + 1, P.Rho * P.Cp) * V
                kf = np.full(N, P.k)
            else:
                M = rho_f(C_mid) * cp_f(C_mid) * V
                kf = arith(k_f(C_mid))
            lo2, di2, up2 = assemble(kf, rf, dr, N, P.h * P.R)
            if not checked:
                check_K(lo2, di2, up2, P.h * P.R, N)
                checked = True
            bT = M / dt * T_old
            bT[-1] += P.R * P.h * T_inf(t1)
            T_next = bsolve(lo2, di2, up2, M / dt, bT, N)

            err = max(np.abs(T_next - T_new).max(),
                      np.abs(C_next - C_new).max())
            T_new, C_new = T_next, C_next
            if err < P.tol:
                break
        else:
            warnings.warn(f'Picard 未收敛：t={t1:.3f}s，残差={err:.3e}，'
                          f'已达上限 maxit={P.maxit}')
        itmax = max(itmax, it + 1)

        # 离散守恒残差：验证组装与线性求解，不构成精度证据（PROJECT.md §12 第9条）
        mmax = max(mmax, abs((V * (C_new - C_old)).sum() / dt
                             - P.hm * P.R * (C_inf(t1) - C_new[-1])))
        emax = max(emax, abs((M * (T_new - T_old)).sum() / dt
                             - P.R * P.h * (T_inf(t1) - T_new[-1])))

        T, C = T_new, C_new
        if out_times is None or round(t1, 6) in out_times:
            out[t1] = (T.copy(), C.copy())
        if verbose and n % log_every == 0:
            print(f'    t={t1 / 3600:.1f}h  T中心={T[0]:.3f}  C最大={C.max():.4f}')

    return out, {'r': r, 'V': V, 'itmax': itmax, 'emax': emax, 'mmax': mmax}


if __name__ == '__main__':
    T_inf, C_inf = load_env()

    print('=' * 72)
    print('检验1  退化回归：冻结为问题1 的常物性')
    out, diag = solve(1800.0, T_inf, C_inf, N=400, dt=0.5,
                      q1_mode=True, out_times={1800.0})
    T, C = out[1800.0]
    REF_T, REF_C = 33.5753, 1.5102          # 问题1 正式结果（outputs/q1/data）
    print(f'  t=1800s  中心T={T[0]:.4f}(问1 {REF_T})   表面C={C[-1]:.4f}(问1 {REF_C})')
    print(f'  偏差：T {abs(T[0] - REF_T):.4f} °C，C {abs(C[-1] - REF_C):.4f} kg/kg')
    print('  说明：问题1 用 Crank-Nicolson，本求解器用后向欧拉，残余偏差为时间格式差异。')

    print('\n检验2  网格与时间步收敛 (t=3h)')
    for Nx, dtx in [(100, 2.0), (200, 1.0), (400, 0.5)]:
        o, _ = solve(P.t_end, T_inf, C_inf, N=Nx, dt=dtx, out_times={P.t_end})
        T, C = o[P.t_end]
        print(f'  N={Nx:4d} dt={dtx:4.1f}  中心T={T[0]:.5f}  表面T={T[-1]:.5f}  '
              f'中心C={C[0]:.5f}  表面C={C[-1]:.5f}')

    print('\n求解问题2 ...')
    out, diag = solve(P.t_end, T_inf, C_inf,
                      out_times=np.arange(1, int(P.t_end) + 1, dtype=float),
                      verbose=True)
    print(f"  Picard 最大迭代 = {diag['itmax']}")
    print(f"  能量守恒绝对残差 = {diag['emax']:.2e}")
    print(f"  质量守恒绝对残差 = {diag['mmax']:.2e}")

    idx = np.arange(0, P.N + 1, P.N // 20)
    r_out = diag['r'][idx] * 100
    times = np.arange(1, int(P.t_end) + 1)
    T_all = np.array([out[float(t)][0][idx] for t in times])
    C_all = np.array([out[float(t)][1][idx] for t in times])

    # 物理自检：极值原理 + 表面单调失水 + 表面不湿于中心
    T_env_max = float(np.max(T_inf(times.astype(float))))
    assert (T_all >= min(P.T0, float(np.min(T_inf(times.astype(float))))) - 1e-6).all(), \
        '温度低于极值原理下界'
    assert (T_all <= max(P.T0, T_env_max) + 1e-6).all(), '温度高于极值原理上界'
    assert np.all(np.diff(C_all[:, -1]) < 1e-9), '表面含水率非单调下降'
    assert (C_all[:, -1] < C_all[:, 0]).all(), '表面应比中心干'
    print('  物理自检通过')

    cols = [0, 5, 10, 15, 20]
    print('\n表3  3小时内药材的温度 (°C)')
    print('时间/h |' + ''.join(f'{r_out[c]:>10.1f}' for c in cols))
    for th in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        k = int(th * 3600) - 1
        print(f'{th:6.1f} |' + ''.join(f'{T_all[k, c]:10.4f}' for c in cols))

    print('\n表4  3小时内药材的水分浓度 (kg/kg)')
    print('时间/h |' + ''.join(f'{r_out[c]:>10.1f}' for c in cols))
    for th in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        k = int(th * 3600) - 1
        print(f'{th:6.1f} |' + ''.join(f'{C_all[k, c]:10.4f}' for c in cols))

    os.makedirs(DEFAULT_OUTPUT_DATA, exist_ok=True)
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q2_T.npy'), T_all)
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q2_C.npy'), C_all)
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q2_r.npy'), r_out)
    print(f'\n中间结果已保存到 {DEFAULT_OUTPUT_DATA}')
    print('待办：result2.xlsx 仍需按问题1 的规格从官方模板导出并回读验收。')
