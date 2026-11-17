# Участие в разработке Bell

Bell — пакет Julia для симплектического моделирования пучка в ускорительной решётке.
Ниже — процесс для людей. Правила для AI-агентов — в [`AGENTS.md`](AGENTS.md).

## GitHub Flow

Цепочка одна: **Issue → branch → PR → CI green → squash-merge → main**.

1. Откройте Issue (bug, feature или physics-task) либо возьмите уже открытый.
2. Ветка от актуального `main`.
3. Один pull request закрывает этот Issue.
4. CI зелёный: тесты и сборка документации.
5. В `main` изменения попадают только через squash-merge.

## Branch naming

Префикс ветки — тип задачи:

| префикс | когда |
| --- | --- |
| `feat/` | новая возможность |
| `fix/` | исправление |
| `refactor/` | то же поведение, другая структура кода |
| `test/` | только тесты |
| `docs/` | только документация |
| `chore/` | инфраструктура, CI, шаблоны |

Пример: `feat/thick-rf-map`, `fix/drift-hamiltonian`.

## Conventional Commits

Заголовок коммита:

```text
<type>: <краткое описание>
```

Типы те же, что у веток: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`.
В теле коммита пишите зачем, а не список файлов.

Squash-merge оставляет в `main` один коммит на PR, поэтому заголовок PR —
тоже conventional commit.

## Squash-merge

- В `main` сливается только squash-merge. Merge-коммиты и rebase-merge не используем.
- Один PR = одна задача = один Issue.
- Несвязанные правки раскладывайте по разным PR.

## Tests

Перед коммитом и перед PR:

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

CI гоняет то же самое на Julia `1.10` и текущей стабильной Julia `1`,
на `ubuntu-latest` и `windows-latest`.

## Versioning

Не повышайте `version` в `Project.toml` в обычном PR.
Версию меняют только в релизном PR, когда публикуется новая версия пакета.

## CHANGELOG

Журнал ведётся по [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Канонический файл — [`CHANGELOG.md`](CHANGELOG.md) в корне репозитория.

В пользовательском PR допишите заметку в секцию `[Unreleased]`
(`Added`, `Changed`, `Fixed` или `Removed`).
Правку, которую не видит пользователь пакета (только CI или шаблоны),
в секцию можно не писать.

Подробные заметки до версии 0.17 лежат в архиве [`docs/HISTORY.md`](docs/HISTORY.md).
Архив не удаляем и не переписываем.

## AI agents

Grok, Cursor и другие агенты следуют [`AGENTS.md`](AGENTS.md):
один сеанс — один Issue, одна ветка, один PR.
