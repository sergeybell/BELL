# Changelog

Все заметные изменения Bell записываются в этот файл.

Формат — [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
версии — [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Подробные исторические заметки (v0.1–v0.17) остаются в
[`docs/HISTORY.md`](docs/HISTORY.md). Это архив, его не переносят сюда покоммитно.
Дальше канонический журнал — этот файл.

Краткий указатель тех же версий по-прежнему лежит в
[`docs/CHANGELOG.md`](docs/CHANGELOG.md).

## [Unreleased]

### Added

- `ParticleState` (канонический 6-вектор и слот спина, по умолчанию пустой) и `Beam` (набор таких частиц). Уравнение Томаса–Баргмана–Мишеля–Телегди не считается: спин только хранится.
- `AbstractMap`, обёртки `LinearMap` и `DAMap` (`linear_part`, вызов как функция, `compose` для пары обёрток). `LieFactor` — тоже `AbstractMap`.
- Трекеры `FieldTracker` (`exact_track_turn`), `LieTracker` (`track` по факторам `build_factors`) и `MapTracker` (однооборотная функция), плюс `track!` для `Beam` и для изменяемого 6-вектора. Старые функции трекинга остаются на месте.
- Тесты разложены по файлам `test/test_*.jl`: round-trip структуры и координат, регрессия оптики `magnetic_2` по `benchmarks/twiss/bell/results.csv`, проверка Aqua.
- `test/test_physics_consistency.jl`: `element_field` совпадает с `S ∇H` того же гамильтониана на DA порядков 5 и 6 и с численным градиентом на `Float64`.
- README: установка через `Pkg.develop`, быстрый старт, запуск тестов, ссылки на `CONTRIBUTING.md`, `CHANGELOG.md` и сборку документации.
- В начале `docs/HISTORY.md` указано, что это архив до v0.17, а новые записи ведутся в этом файле.

### Changed

- CI — одна задача: `windows-latest`, Julia 1. Матрица ОС и версий, сбор покрытия, загрузка в Codecov и сборка документации убраны. Задача запускает обычный `Pkg.test()`.
- `Pkg.test()` по умолчанию гоняет быстрый набор. Регрессия оптики, DA высоких порядков, Aqua и тесты оптимизаторов запускаются только при `BELL_FULL_TESTS=1`.
- `LieFactor` является подтипом `AbstractMap`. Вызов `F(z)` делает то же, что `apply_factor(F, z)`.
- HTML-документация собирается из одного дерева `docs/src`. Страницы сайта перенесены из `docs/site`. Markdown руководств PDF лежит в `docs/manual`.
- Страница архитектуры описывает слои `math → core → io → physics → maps/{lie,da} → analysis → tracking interfaces → effects` (каталога `effects` ещё нет). Оба движка (Lie на TaylorSeries и DA на DiffAlgebra) остаются.
- Публичные имена Bell перечислены явным `export` вместо цикла `@eval export`. Набор имён тот же; добавлена `clear_caches!()`.
- Кэш карт DA (`_DACACHE`) включает элемент, `p0`, порядок, сам метод интегрирования, `tol` и `nsteps`. `clear_caches!()` очищает кэши карт Lie, DA и краевого поля.
- Графики (`plot_twiss`, `plot_dispersion2`, `plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp`, `gif_elements`) перенесены в расширение `BellPlotsExt`. Plots — слабая зависимость: ядро Bell её не требует, методы появляются после `using Plots`.
- `hamiltonian(e, ref, z)` считает тот же s-гамильтониан на числах, `TaylorN` и DA. `hamiltonian(e, ref)` по-прежнему возвращает полином Lie-пути.
- `src/` разложен по слоям `math`, `core`, `io`, `physics`, `maps`, `analysis`. Имена модулей и публичный API не изменились.
- Цепочки `isa` у `element_lie_map`, `taylor_map_element`, `element_da_map`, `exact_track_element`, `linear_matrix`, `element_rhs`, `Base.show` и `_card_line` заменены методами. `magnetic_lattice` отбирает RF через `strip_rf`. Поведение то же; `element_rhs` для `Marker` и `RFCavity` бросает `ArgumentError`.
- `track(factors, z; turns)` больше не принимает игнорируемые `ref` и `path_length`. Симплектичность матрицы, карты `TaylorN`, карты DA и шага трекинга — одна функция `symplectic_error` (случай DA считает `da_symplectic_error`).

### Fixed

### Removed

- Из `DiffAlgebra` удалены неиспользуемые `poisson_bracket`, `lie_series`, `lie_flow`, `vector_field_flow` и `change_order`. Оставлены `canonical_pairs`, `hamiltonian_vector_field` и `da_symplectic_error`.

## [0.17.0] - 2026-11-02

### Added

- s-гамильтониан `RFCavity`: тонкая карта (`L = 0`) совпадает с `rf_kick`; толстая карта использует длину `L` (`exp(-L:H:)` на Lie-пути, drift/kick/drift на field/DA).
- Функции `rf_potential`, `rf_wavenumber` и `magnetic_lattice`.

### Changed

- Раскладка исходников: `src/lie`, `src/da`, `src/common`, `src/longitudinal_models`.
- Карты элемента (`element_lie_map`, tracking, fringe, DA) больше не содержат отдельной ветки для RF.
