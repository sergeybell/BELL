# ВЧ-резонатор и ускорение (v0.17)

```@meta
CurrentModule = Bell
```

## 1. Что реализовано

`RFCavity(name, V, h, phi_s)` — тонкий (L = 0), `RFCavity(name, V, h, phi_s, L)` — толстый.
Резонатор имеет s-гамильтониан, как все элементы (`src/lie/hamiltonians.jl`):

```math
H = \mathrm{energy\_term}(\delta) - p_s + \frac{G}{L},\qquad
G(l,\delta) = \frac{ZV}{p_0c\,k}\cos\!\Big(\varphi_s + \frac{k\,l}{\beta_0\,v(\delta)/v_0}\Big),\qquad
k = \frac{h\omega_0}{c} = \frac{2\pi h\beta_0}{C}.
```

`G = rf_potential(cav, ref, C)`, `k = rf_wavenumber(cav, ref, C)`; V в МВ, p0c в МэВ, `C = path_length(lattice)`.
Модель — однородный зазор с фактором пролётного времени T = 1 (V — эффективное напряжение, V = E_s L T).

**Вывод и знаки.** В паре ``(q,P)=(-c\Delta t,\ E/(p_0c))`` сила зазора
``dP/ds=-\partial H/\partial q = ZV\sin(\varphi_s+kq)/(p_0cL)``. Пара Bell ``(l,\delta)`` получается точечным
преобразованием ``\delta=\delta(P)``, ``l=q\,dP/d\delta=q\beta=-v\Delta t``; гамильтониан не меняется,
``q=l/(\beta_0 v/v_0)``, а ``kq=2\pi h l/(C v/v_0)`` — ровно фаза `rf_kick` и `canonical_to_laboratory`
(``\varphi>0`` — частица пришла раньше). Знак проверен тестами: тонкий предел совпадает с `rf_kick`,
``Q_s`` из однооборотной матрицы совпадает с формулой `Smooth.synchrotron_tune` (устойчиво при
``\eta_0\cos\varphi_s>0``).

**Отображения.** Тонкий: ``\exp(-{:}G{:})`` (в ``(q,P)`` ``G`` зависит только от ``q`` ⇒ поток точный:
``\Delta t`` фиксировано, ``E\mathrel{+}=ZV\sin(\varphi_s+\varphi)`` — это `rf_kick`). Толстый:
``\exp(-L{:}H{:})`` (Lie), в полевом/DA-пути — расщепление Странга «дрейф L/2n – толчок V/n – дрейф L/2n».
`Tracking.build_factors / track / ring_lie_factor / taylor_map_element` и `DATracking.element_da_map(s)`
больше не имеют отдельных веток для RF. Для магнитных величин (хроматичность, α, Twiss) —
`magnetic_lattice(lattice)` (тонкий RF удаляется, толстый заменяется дрейфом).

Проверки (test/runtests.jl, «RF Hamiltonian»): Taylor-карта ``\exp(-{:}G{:})`` порядка 3 vs `rf_kick`
на TaylorN — 2.2e-15; толстый vs «D(L/2)·толчок·D(L/2)» — 7.6e-7, 7.6e-9, 7.9e-11 при L = 0.1, 0.01, 0.001
(O(L²) по коэффициентам); `symplectic_error` 2.6e-23 (тонкий), ≤4.4e-16 (толстый);
``Q_s`` карта / Smooth − 1 = 1.5e-5; 10 оборотов: факторы элементов vs однооборотный фактор 6.5e-12,
Lie vs точное полевое трекирование 5.4e-7.

## 2. Как это сделано в COSY INFINITY (идеология, по cosy.fox, процедура `RF`)

- Резонатор — **тонкий** элемент: `RF V N W PHI D` (V — полином по x, y порядка N, т.е. напряжение может
  зависеть от поперечного положения; W — частота; PHI — фаза в градусах).
- Время прибытия берётся из продольной координаты карты, ``t = l\cdot(-p_0/E_0/c)``; толчок
  ``\Delta E = Z\,V(x,y)\sin(2\pi W t + \varphi)`` (с учётом отличия массы/заряда частицы от опорной).
- Опорная энергия **сразу меняется** на набор синхронной частицы: `RP E0+VV0 M0 Z0`. После этого
  поперечные импульсы масштабируются на ``p_{0,\rm old}/p_{0,\rm new}`` (адиабатическое затухание),
  продольная координата пересчитывается из сохранённого времени с новыми ``E_0,p_0``,
  а энергетическое отклонение перенормируется: ``\delta_K' = (\delta_K + (\Delta E-\Delta E_0)/E_0)/(1+\Delta E_0/E_0)``.
- То есть COSY не строит гамильтониан резонатора: это явное (DA-)отображение «толчок + смена опоры»,
  применяемое к карте. Его якобиан в нормированных переменных несимплектичен (det = (p_old/p_new)^k) —
  это и есть затухание.

**Чем Bell отличается.**
1. Резонатор — гамильтонов элемент; тонкий толчок получается как предел, а не постулируется; есть
   толстый резонатор с тем же H; все отображения (Lie, Taylor, DA) строятся одним путём.
2. Накопитель (``\sin\varphi_s=0``): опора не меняется, отображение симплектично и не зависит от оборота —
   можно строить однооборотную карту, нормальную форму 6D, ``Q_s``.
3. Ускорение: смена опоры отделена от физики резонатора и описана как каноническое преобразование
   (см. §3); пока используется дискретный `renormalize` (как COSY и Xsuite), но с точным пересчётом через
   ``\Delta t`` и ``v(\delta)``, в переменных ``\delta = \Delta p/p_0``, а не ``\delta_K``.
4. Bell не поддерживает пока ``V(x,y)`` и разный заряд/массу частиц (в COSY есть) — оставлено.

## 3. Ускорение в гамильтониане — TODO (вывод)

При ``\sin\varphi_s\ne0`` у ``H`` есть линейный член по ``l``: ``-\partial H_{rf}/\partial l|_0 = ZV\sin\varphi_s/(\beta_0p_0cL)``.
Вводим опорный импульс ``p_0(s)`` с ``dE_0/ds = ZV\sin\varphi_s/L`` (``dE_0 = v_0dp_0``). Производящая функция
``F_2 = p_0(s)[x p_x + y p_y + l(1+\delta)]`` даёт добавку

```math
\Delta H = \frac{p_0'}{p_0}\big[x p_x + y p_y + l(1+\delta)\big],\qquad
-\frac{\partial \Delta H}{\partial l}\Big|_{\delta=0} = -\frac{p_0'}{p_0} = -\frac{ZV\sin\varphi_s}{\beta_0p_0cL},
```

т.е. линейный член сокращается, ``l=\delta=0`` снова решение, а ``x p_x`` даёт адиабатическое затухание.
При этом ``\beta_0(s),\gamma_0(s)`` в `energy_term` и ``k(s)`` зависят от ``s`` — гамильтониан явно зависит от ``s``,
однооборотная карта зависит от оборота; магнитные поля при этом считаются растущими с ``p_0``
(нормированные k1, k2, h постоянны). Реализация требует `RefParticle`, зависящего от s внутри элемента
(или шаговой смены опоры на срезах) — не сделано в v0.17. Сейчас работает путь `Acceleration.accelerate`
(`src/longitudinal_models/`): `rf_kick` + `renormalize` после каждого резонатора (дискретная форма той же
замены, Xsuite §2.21: ``\delta_{new}=(1+\delta_{old})p_{0,old}/p_{0,new}-1``, ``\zeta_{new}=\zeta_{old}\beta_{0,new}/\beta_{0,old}``).

## 4. Литература

- A. Wolski, Lecture 4 «RF cavities», Cockcroft Institute (2014), ур. 31–37 — H моды TM010, усреднение по
  ``L=\pi/k``, фактор T, тонкий предел ``\delta_W \mathrel{+}= (qV/P_0c)\sin(\varphi_0-kz)``.
  https://www.cockcroft.ac.uk/wp-content/uploads/2014/12/wolski-4.pdf
  (у Вольского ``\delta_W = E/(p_0c) - 1/\beta_0`` — энергия; ``z = q``; у Bell ``1+\delta=\sqrt{(\delta_W+1/\beta_0)^2-1/(\beta_0\gamma_0)^2}``, ``l = \beta z``.)
- A. Wolski, Lecture 7 «Longitudinal dynamics», Cockcroft Institute (2014).
- Xsuite physics manual, §2.14 (Cavity), §2.21 (изменение опорной энергии). https://xsuite.readthedocs.io
- F. Tecker, «Longitudinal beam dynamics», CAS 2020, arXiv:2011.02932.
- J. Le Duff, «Longitudinal beam dynamics in circular accelerators», CAS, CERN 94-01.
- COSY INFINITY 10.2, `cosy.fox`, процедура `RF` (прочитана только идеология, код не копировался).
