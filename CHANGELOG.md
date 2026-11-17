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

### Changed

### Fixed

### Removed

## [0.17.0] - 2026-11-02

### Added

- s-гамильтониан `RFCavity`: тонкая карта (`L = 0`) совпадает с `rf_kick`; толстая карта использует длину `L` (`exp(-L:H:)` на Lie-пути, drift/kick/drift на field/DA).
- Функции `rf_potential`, `rf_wavenumber` и `magnetic_lattice`.

### Changed

- Раскладка исходников: `src/lie`, `src/da`, `src/common`, `src/longitudinal_models`.
- Карты элемента (`element_lie_map`, tracking, fringe, DA) больше не содержат отдельной ветки для RF.
