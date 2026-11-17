# Bell

**BELL — Beam Evolution in Accelerator Lattice**

**Автор:** С. Д. Колокольчиков (Kolokolchikov Sergey Dmitrievich, `sergey.bell13@gmail.com`)

Bell — пакет Julia для симплектического моделирования пучка в ускорительной решётке.
Гамильтониан элемента становится преобразованием Ли; произведение этих преобразований — объект одного оборота.
Из него читаются оптика и величины старших порядков. Точное интегрирование поля — независимая проверка.
Оба движка остаются в пакете: путь Ли на TaylorSeries и путь DA на DiffAlgebra.

## Установка

В ту среду Julia, где нужен Bell (`path` — каталог клона с `Project.toml` пакета):

```julia
using Pkg
Pkg.develop(path = raw"C:\path\to\Bell")
```

После этого `using Bell` берёт этот клон: правки в дереве видны без переустановки.

Зависимости самого пакета из корня репозитория:

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Примеры. В `examples/Project.toml` путь к Bell задан таблицей `[sources]` (Julia 1.11 и новее):

```text
julia --project=examples -e "using Pkg; Pkg.instantiate()"
```

На Julia 1.10 ту же связь задают явно, из корня репозитория:

```text
julia --project=examples -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
```

## Быстрый старт

Из корня репозитория:

```julia
using Bell

lat = load_structure("examples/lattices/fodo.bell")
R = ring_matrix(lat)
opt = periodic_optics(R, lat.ref, path_length(lat))
opt.x.beta, opt.y.beta
```

Линейная оптика решётки `magnetic_2`:

```text
julia --project=examples examples/10_linear_optics/run.jl
```

## Тесты

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

## Документация и процесс

* Сайт Documenter, одно дерево [`docs/src/index.md`](docs/src/index.md). Сборка совпадает с задачей CI:

  ```text
  julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
  julia --project=docs docs/make.jl
  ```

  На Julia 1.11 и новее достаточно `Pkg.instantiate()` в проекте `docs`: путь к Bell уже записан в `docs/Project.toml`.
* Руководство по физике (PDF): Markdown в [`docs/manual/index.md`](docs/manual/index.md), сборка `julia --project=docs docs/make_pdf.jl`. Готовый PDF — `docs/pdf/Bell_Beam_Physics_Manual.pdf`.
* Журнал изменений: [`CHANGELOG.md`](CHANGELOG.md). Заметки до v0.17 — архив [`docs/HISTORY.md`](docs/HISTORY.md). Краткий указатель тех же версий — [`docs/CHANGELOG.md`](docs/CHANGELOG.md).
* Как предлагать изменения: [`CONTRIBUTING.md`](CONTRIBUTING.md). Правила для агентов: [`AGENTS.md`](AGENTS.md).
