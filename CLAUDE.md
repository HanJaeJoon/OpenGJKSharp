# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Commit messages

- This is an open-source project: **write all commit messages in English.**
- Format: `<type>: <subject>` (single-service convention; no `[service]` prefix).
  - Types: `feat`, `fix`, `refactor`, `style`, `docs`, `test`, `chore`, `perf`, `ci`, `hotfix`
- Subject: imperative mood, 50 characters or less, no trailing period.
- Body (optional): explain what and why, in English.

## Project notes

- Library source: `src/OpenGJKSharp`, tests: `tests/OpenGJKSharp.Tests` (xunit).
- Multi-targeted package: `netstandard2.0;netstandard2.1;net8.0;net10.0`.
  Keep new code compatible with netstandard2.0 (no `Math.Clamp`, `float.IsFinite`, etc. -- use the internal helpers).
- NuGet publishing runs only on `v*` tag push and verifies the tag matches `<Version>` in the csproj.
- Never run `dotnet pack` on its own: `GeneratePackageOnBuild=true` makes it skip the build and
  package a stale `bin` (or fail with NU5026 when `bin` is clean). Use `dotnet build -c Release`,
  or `dotnet build` followed by `dotnet pack --no-build` -- the order the CI workflow uses.
- At each release, bump `<PackageValidationBaselineVersion>` in the csproj to the latest published version.
- Keep `CHANGELOG.md` (Keep a Changelog format) up to date with user-visible changes.
