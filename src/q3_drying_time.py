"""
2026 高教社杯 A题 —— 问题3 求解器
在问题2 的耦合模型上加终止判据：各处 C < 0.15 kg/kg

相对问题2 的三处变更：
  1) 环境边界需外推到数十小时（附件1 只到 4h），取 t>=9000s 的平台均值
  2) 变步长：早期 1s（预热陡变），后期 60s（干燥缓慢），均整除60以对齐输出
  3) 网格加密到 N=800：终态表面存在 D 坍塌 256 倍的"硬壳层"，粗网格会
     严重高估干燥时间（N=100 给 71.4h，N=800 给 57.5h）

入库说明（2026-09-12）：本文件由项目外的 solve_Q3.py 并入，只改动文件路径
以适配项目结构（附件从 data/official 读取，.npy 写入 outputs/q3/data），
离散格式、步长调度与终止判据一行未动。

待办：
  * docs/q3/drying_time_independent_review.md 指出界面调和平均在 D 退化区
    会低估表面通量（方法 22 的“干不完”为数值假象），改动界面系数时需同步
    更新 docs/q3/drying_time_method.md。
"""
import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from q2_coupled import (P, rho_f, cp_f, k_f, D_f, make_grid, assemble,
                        bsolve, harm, arith, check_K)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ENV_PATH = os.path.join(PROJECT_ROOT, 'data', 'official', 'A题',
                                '附件', '附件1.xlsx')
DEFAULT_OUTPUT_DATA = os.path.join(PROJECT_ROOT, 'outputs', 'q3', 'data')
C_TARGET = 0.15          # kg/kg 干燥终止判据
OUT_DT   = 60.0          # s     result3 输出间隔


def load_env_ext(path=None, plateau_t=9000.0):
    """分段线性插值 + 恒温干燥段常值延拓（取平台均值，避开末点噪声）"""
    if path is None:
        path = DEFAULT_ENV_PATH
    d = pd.read_excel(path).values.astype(float)
    t, T, C = d[:, 0], d[:, 1], d[:, 2]
    m = t >= plateau_t
    Te, Ce = T[m].mean(), C[m].mean()
    return (interp1d(t, T, 'linear', bounds_error=False, fill_value=(T[0], Te)),
            interp1d(t, C, 'linear', bounds_error=False, fill_value=(C[0], Ce)),
            Te, Ce)


def dt_schedule(t):
    """变步长。取值必须整除 OUT_DT=60，保证输出时刻精确落在步点上。"""
    if t < 600:     return 1.0      # 预热初期：边界阶跃、梯度最陡
    if t < 3600:    return 5.0      # 预热中期
    if t < 14400:   return 20.0     # 预热结束、进入恒温段
    return 60.0                     # 恒温干燥期：R^2/D ~ 140h，步长可放大


def solve_Q3(N=800, tMax=150*3600.0, verbose=True):
    Tinf, Cinf, Te, Ce = load_env_ext()
    dr, r, V, rf = make_grid(N, P.R)

    T = np.full(N + 1, P.T0)
    C = np.full(N + 1, P.C0)
    t = 0.0
    idx = np.arange(0, N + 1, N // 20)          # 0,0.1,...,2.0 cm
    rec_t, rec_C = [], []                       # 每 60s 的输出
    t_end, itmax, m_res = None, 0, 0.0
    checked = False

    while t < tMax:
        dt = dt_schedule(t)
        t0, t1 = t, t + dt
        To, Co = T.copy(), C.copy()
        Tn, Cn = T.copy(), C.copy()

        for it in range(P.maxit):
            # (a) 传质：D 在时间中点取值
            Dm = D_f(0.5 * (Co + Cn), 0.5 * (To + Tn))
            lo, di, up = assemble(harm(Dm), rf, dr, N, P.hm * P.R)
            if not checked:
                check_K(lo, di, up, P.hm * P.R, N)
            bC = V / dt * Co
            bC[-1] += P.hm * P.R * Cinf(t1)
            C2 = bsolve(lo, di, up, V / dt, bC, N)

            # (b) 传热：物性用刚更新的 C
            Cm = 0.5 * (Co + C2)
            M = rho_f(Cm) * cp_f(Cm) * V
            lo2, di2, up2 = assemble(arith(k_f(Cm)), rf, dr, N, P.h * P.R)
            if not checked:
                check_K(lo2, di2, up2, P.h * P.R, N); checked = True
            bT = M / dt * To
            bT[-1] += P.R * P.h * Tinf(t1)
            T2 = bsolve(lo2, di2, up2, M / dt, bT, N)

            err = max(np.abs(T2 - Tn).max(), np.abs(C2 - Cn).max())
            Tn, Cn = T2, C2
            if err < P.tol:
                break
        itmax = max(itmax, it + 1)

        m_res = max(m_res, abs((V * (Cn - Co)).sum() / dt
                               - P.hm * P.R * (Cinf(t1) - Cn[-1])))
        T, C = Tn, Cn
        t = t1

        if abs(t % OUT_DT) < 1e-9:
            rec_t.append(t); rec_C.append(C[idx].copy())
        if t_end is None and C.max() < C_TARGET:
            t_end = t
            if abs(t % OUT_DT) > 1e-9:          # 补齐到整 60s
                rec_t.append(t); rec_C.append(C[idx].copy())
            break
        if verbose and abs(t % (6 * 3600)) < 1e-9:
            print(f'    t={t/3600:5.1f}h  C_max={C.max():.5f}  '
                  f'C表面={C[-1]:.5f}  T中心={T[0]:.3f}')

    return dict(t_end=t_end, t=np.array(rec_t), C=np.array(rec_C),
                r_out=r[idx] * 100, itmax=itmax, m_res=m_res,
                C_fin=C, T_fin=T, Te=Te, Ce=Ce, r=r)


if __name__ == '__main__':
    print('=' * 70)
    print('问题3：求 C<0.15 kg/kg 的干燥时间')
    print()
    Tinf, Cinf, Te, Ce = load_env_ext()
    print(f'恒温干燥段延拓：T∞={Te:.4f} °C，C∞={Ce:.5f} kg/kg（t≥9000s 平台均值）')
    print(f'终止时刻的传质驱动力 C_target − C∞ = {C_TARGET - Ce:.5f}')
    print()

    print('网格收敛性（干燥时间 / h）：')
    conv = {}
    for Nx in [200, 400, 800]:
        rr = solve_Q3(N=Nx, verbose=False)
        conv[Nx] = rr['t_end'] / 3600
        print(f'  N={Nx:4d}  {conv[Nx]:.3f} h')
    d1, d2 = conv[200] - conv[400], conv[400] - conv[800]
    print(f'  Richardson 外推(p=2): {conv[800]-d2/3:.3f} h   '
          f'(N800 相对外推值偏差 {abs(d2/3):.3f} h)')
    print()

    print('正式求解 (N=800) ...')
    res = solve_Q3(N=800)
    te = res['t_end']
    print()
    print(f"  干燥结束时间 = {te:.0f} s = {te/3600:.2f} h = {te/86400:.3f} 天")
    print(f"  Picard 最大迭代 = {res['itmax']}")
    print(f"  质量守恒绝对残差 = {res['m_res']:.2e}")
    print(f"  结束时刻 C 范围 [{res['C_fin'].min():.5f}, {res['C_fin'].max():.5f}]")
    print(f"  结束时刻 T 范围 [{res['T_fin'].min():.4f}, {res['T_fin'].max():.4f}]")

    # 物理自检
    C = res['C']
    assert np.all(np.diff(C[:, 0]) < 1e-9),  '中心含水率非单调下降'
    assert np.all(C[:, -1] <= C[:, 0] + 1e-9), '表面应始终不湿于中心'
    assert C[-1].max() < C_TARGET,           '终止时刻未全部达标'
    assert C[-2].max() >= C_TARGET if len(C) > 1 else True, '终止时刻不是首次达标'
    print('  物理自检通过')

    os.makedirs(DEFAULT_OUTPUT_DATA, exist_ok=True)
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q3_t.npy'), res['t'])
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q3_C.npy'), res['C'])
    np.save(os.path.join(DEFAULT_OUTPUT_DATA, 'Q3_r.npy'), res['r_out'])

    # ---- 表5：每 6h + 结束时刻，每 0.5cm ----
    cols = [0, 5, 10, 15, 20]
    print('\n表5  药材烘干过程的水分浓度 (kg/kg)')
    print('时间/h |' + ''.join(f"{res['r_out'][c]:>10.1f}" for c in cols))
    marks = list(np.arange(6, te / 3600, 6)) + [te / 3600]
    for th in marks:
        k = int(round(th * 3600 / OUT_DT)) - 1
        k = min(k, len(res['t']) - 1)
        lab = f'{th:6.2f}' if abs(th - te / 3600) < 1e-9 else f'{th:6.0f}'
        print(lab + ' |' + ''.join(f"{res['C'][k,c]:10.4f}" for c in cols))
