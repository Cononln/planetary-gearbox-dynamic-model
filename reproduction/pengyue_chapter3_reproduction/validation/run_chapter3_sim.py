# -*- coding: utf-8 -*-
"""
彭悦第3章全部振动信号图复现入口。

工况 (21): healthy; sun q3g45; sun q=2/4/6/8@g45; sun g=15/30/45/60/75@q4;
           planet q3g45(P1); planet q=2/4/6/8@g45; planet g=15/30/45/60/75@q4
两个误差版本:
  e20: 论文文字模型, e_av=20μm 全量进入啮合 (图3-9时域含1fm大幅分量)
  e0 : 论文图形自洽版本 (论文图3-10的1fm幅度与e_av=20μm响应矛盾约45倍,
        判定其绘图信号不含1fm误差准静态分量), 用于与论文图对比
输出: figures/fig3_*.png, results/chapter3_sim/*.csv
"""
import sys, os, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance
from pengyue_chapter3_reproduction.dynamics.model18 import Model18, rk4_run

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = os.path.join(os.path.dirname(__file__), '..')
FIG = os.path.join(ROOT, 'figures')
RES = os.path.join(ROOT, 'results', 'chapter3_sim')
os.makedirs(FIG, exist_ok=True); os.makedirs(RES, exist_ok=True)

EPS_SP = P.contact_ratio('sp')
EPS_PR = P.contact_ratio('pr')
ETA = np.linspace(0, 1, 501)
A_SP = P.zs * P.phi_p / (2 * np.pi)
A_RP = -P.zr * P.phi_p / (2 * np.pi)
SUN_FAULT_TOOTH = 15
PLANET_FAULT_TOOTH = 15
PLANET_FAULT_IDX = 0
DT = (1.0 / P.f_m) / 64
T_REC = 16 * P.zs / P.f_m          # 1.4242 s 记录段
T_TRAN = 2.0                        # 瞬态舍弃(最慢支撑模态~20Hz, ζ=0.03)

_tbl_cache = {}

def sp_tables(q_mm=None, g_deg=None, member='sun'):
    key = ('sp', member, q_mm, g_deg)
    if key in _tbl_cache:
        return _tbl_cache[key]
    crack = dict(q=q_mm * 1e-3, gamma=np.deg2rad(g_deg)) if q_mm else None
    ch = external_pair_compliance(ETA, n_int=300)
    if crack:
        if member == 'sun':
            cc = external_pair_compliance(ETA, crack_sun=crack, n_int=300)
        else:
            cc = external_pair_compliance(ETA, crack_planet=crack, n_int=300)
    else:
        cc = None
    _tbl_cache[key] = (ch, cc)
    return _tbl_cache[key]

def pr_tables(q_mm=None, g_deg=None):
    key = ('pr', q_mm, g_deg)
    if key in _tbl_cache:
        return _tbl_cache[key]
    crack = dict(q=q_mm * 1e-3, gamma=np.deg2rad(g_deg)) if q_mm else None
    ch = internal_pair_compliance(ETA, n_int=300)
    cc = internal_pair_compliance(ETA, crack_planet=crack, n_int=300) if crack else None
    _tbl_cache[key] = (ch, cc)
    return _tbl_cache[key]

def pair_curve(comp):
    return 1.0 / comp['c_total']

def build_model(case, e_av):
    if case['type'] == 'sun':
        ch, cc = sp_tables(case['q'], case['g'], 'sun')
        ksp_c = pair_curve(cc) if cc is not None else None
        m = Model18((ETA, pair_curve(ch)), (ETA, pair_curve(pr_tables()[0])),
                    EPS_SP, EPS_PR, A_SP, A_RP,
                    fault_sp_tooth=SUN_FAULT_TOOTH if ksp_c is not None else None,
                    ksp_cracked=ksp_c, friction=False, e_av=e_av)
    elif case['type'] == 'planet':
        ch, cc = sp_tables(case['q'], case['g'], 'planet')
        prh, prc = pr_tables(case['q'], case['g'])
        has = case['q'] is not None
        m = Model18((ETA, pair_curve(ch)), (ETA, pair_curve(prh)),
                    EPS_SP, EPS_PR, A_SP, A_RP,
                    fault_sp_tooth=PLANET_FAULT_TOOTH if has else None,
                    ksp_cracked=pair_curve(cc) if has else None,
                    sp_fault_gear='planet',
                    fault_pr=dict(planet=PLANET_FAULT_IDX, tooth=PLANET_FAULT_TOOTH,
                                  k=pair_curve(prc), offset_cycles=P.zp / 2.0) if has else None,
                    friction=False, e_av=e_av)
    else:
        m = Model18((ETA, pair_curve(sp_tables()[0])), (ETA, pair_curve(pr_tables()[0])),
                    EPS_SP, EPS_PR, A_SP, A_RP, friction=False, e_av=e_av)
    return m

def static_equilibrium(m, t_ref=0.0):
    """牛顿法解静平衡 (e_av 置零, 刚度取 t_ref 时刻), 消除阶跃瞬态。"""
    e_keep = m.e_av
    m.e_av = 0.0
    q = np.zeros(18)
    for _ in range(8):
        g = m.rhs(t_ref, q, np.zeros(18)) * m.M   # 恢复力量纲
        J = np.zeros((18, 18))
        for j in range(18):
            dq = np.zeros(18); dq[j] = 1e-7 * max(1.0, abs(q[j]))
            J[:, j] = (m.rhs(t_ref, q + dq, np.zeros(18)) * m.M - g) / dq[j]
        dq = np.linalg.solve(J, g)
        q = q - dq
        if np.max(np.abs(dq)) < 1e-12:
            break
    m.e_av = e_keep
    return q


def run_case(case, e_ver):
    e_av = 20e-6 if e_ver == 'e20' else 0.0
    m = build_model(case, e_av)
    q0 = static_equilibrium(m)
    T, Q, _ = rk4_run(m, T_TRAN + T_REC, DT, q0=q0, record=True)
    sel = T >= T_TRAN
    t = T[sel] - T_TRAN
    dte = np.stack([np.array([m.dte_sp(t[j] + T_TRAN, Q[j], i) for j in range(t.size)])
                    for i in range(3)])
    dte_pr1 = np.array([m.dte_pr(t[j] + T_TRAN, Q[j], 0) for j in range(t.size)])
    # 动态分量: 扣除准静态啮合变形 F0/k(t) (论文图3-9/3-10结构对应信号, 见报告§2)
    F0 = m.Tin / (3 * m.rbs)
    dte_dyn = np.empty_like(dte)
    for i in range(3):
        kk = np.array([m.k_sp(t[j] + T_TRAN, i) for j in range(t.size)])
        dte_dyn[i] = dte[i] - F0 / kk
    kpr = np.array([m.k_pr(t[j] + T_TRAN, 0) for j in range(t.size)])
    dtepr_dyn = dte_pr1 - F0 / kpr
    return t, dte, dte_pr1, m, dte_dyn, dtepr_dyn

def spectrum(x):
    X = np.fft.rfft(x) / x.size
    amp = 2 * np.abs(X)
    df = 1.0 / (x.size * DT)
    return amp, df, np.arange(amp.size) * df

CASES = [
    dict(name='healthy', type='none', q=None, g=None),
    dict(name='sun_q3_g45', type='sun', q=3, g=45),
] + [dict(name=f'sun_q{q}_g45', type='sun', q=q, g=45) for q in (2, 4, 6, 8)] \
  + [dict(name=f'sun_q4_g{g}', type='sun', q=4, g=g) for g in (15, 30, 60, 75)] \
  + [dict(name='pl_q3_g45', type='planet', q=3, g=45)] \
  + [dict(name=f'pl_q{q}_g45', type='planet', q=q, g=45) for q in (2, 4, 6, 8)] \
  + [dict(name=f'pl_q4_g{g}', type='planet', q=4, g=g) for g in (15, 30, 60, 75)]

def main(runs=('e0',)):
    store = {}
    t0 = time.time()
    for e_ver in runs:
        for case in CASES:
            key = (e_ver, case['name'])
            print(f'--- run {key} ...', flush=True)
            t, dte, dtepr, m, dte_dyn, dtepr_dyn = run_case(case, e_ver)
            store[key] = dict(t=t, dte=dte, dtepr=dtepr,
                              dte_dyn=dte_dyn, dtepr_dyn=dtepr_dyn)
            amp, df, fr = spectrum(dte_dyn[0])
            store[key].update(amp=amp, df=df, fr=fr)
            np.savetxt(os.path.join(RES, f'dte_{e_ver}_{case["name"]}.csv'),
                       np.column_stack([t, dte_dyn[0] * 1e6, dte_dyn[1] * 1e6,
                                        dte[0] * 1e6, dtepr_dyn * 1e6]),
                       delimiter=',',
                       header='t_s,xsp1dyn_um,xsp2dyn_um,xsp3dyn_um,xsp1raw_um,xpr1dyn_um',
                       comments='')
    print(f'All runs done in {time.time()-t0:.0f}s')
    return store

if __name__ == '__main__':
    st = main(runs=('e0', 'e20'))
    np.save(os.path.join(RES, 'store_meta.npy'), np.array([DT, T_REC, T_TRAN]))
    import pickle
    with open(os.path.join(RES, 'sim_store.pkl'), 'wb') as fh:
        pickle.dump(st, fh)
    print('store saved')
