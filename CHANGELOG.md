# phac-nml/profile_dists: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.8] - 2025/06/11

### `Fixed`

- Fixed [issue 36](https://github.com/phac-nml/profile_dists/issues/36) introduced in [1.0.7] where `--column` would cause error. [PR #37](https://github.com/phac-nml/profile_dists/pull/37)

### `Added`

- Added a pytest to test functionality of `--column` and close the [issue 12](https://github.com/phac-nml/profile_dists/issues/12). [PR #37](https://github.com/phac-nml/profile_dists/pull/37)

## [1.0.7] - 2025/06/09

### `Fixed`

- Addressed an issue where certain cases of missing loci were leading to self-distances. [PR #34](https://github.com/phac-nml/profile_dists/pull/34)

## [1.0.6] - 2025/05/21

### Fixed
- Updated many versions in `setup.py`. [PR #30](https://github.com/phac-nml/profile_dists/pull/30)

## [1.0.5] - 2025/05/01

### `Fixed`

- Fix function `update_column_map` to resolve [#25](https://github.com/phac-nml/profile_dists/issues/25). [PR #27](https://github.com/phac-nml/profile_dists/pull/27)

## [1.0.4] - 2025/03/31

### `Fixed`

- Fixed issue with calculating distances in the scenario of different sets of samples in the query and reference files. See [PR #19](https://github.com/phac-nml/profile_dists/pull/19).

[1.0.4]: https://github.com/phac-nml/profile_dists/releases/tag/1.0.4
[1.0.5]: https://github.com/phac-nml/profile_dists/releases/tag/1.0.5
[1.0.6]: https://github.com/phac-nml/profile_dists/releases/tag/1.0.6
[1.0.7]: https://github.com/phac-nml/profile_dists/releases/tag/1.0.7
[1.0.8]: https://github.com/phac-nml/profile_dists/releases/tag/1.0.8
