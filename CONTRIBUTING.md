# Contributing to Bell

Workflow for a change to this repository. What the code is, and the rules for agents, are in [`AGENTS.md`](AGENTS.md).

## Workflow

One issue, one branch, one pull request.

1. Open an issue, or take an open one.
2. Branch from current `main`. The name is `<type>/<short-description>`, for example `fix/drift-hamiltonian`.
3. Open a pull request that closes that issue.
4. The fast tests are green.
5. Squash-merge into `main`.

`main` receives squash-merges only.

## Conventional Commits

Commits follow [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).

```text
<type>[optional scope][optional !]: <description>
```

Types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`, `build`.

A breaking change is marked with `!` after the type (and the scope, if there is one), or with a `BREAKING CHANGE:` footer.

Squash-merge keeps one commit on `main`, and that commit is the pull request title. The title is a conventional commit. The body says why.

## Tests

Fast tests, from the repository root:

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

Full tests set `BELL_FULL_TESTS=1` and run the same command. In PowerShell:

```text
$env:BELL_FULL_TESTS = "1"
julia --project=. -e "using Pkg; Pkg.test()"
```

CI on a pull request runs the fast suite: Julia 1 on `windows-latest`.

## Changelog and version

A user-visible change goes in [`CHANGELOG.md`](CHANGELOG.md) under `[Unreleased]`, in `Added`, `Changed`, `Fixed`, or `Removed`. The entry cites the pull request number. The file follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

A change a package user does not see, such as CI or a template, needs no entry.

`version` in `Project.toml` changes only in a release pull request. That pull request renames `[Unreleased]` to the new version and date.
