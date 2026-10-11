## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* NOTE: "File which contains pragma(s) suppressing diagnostics:
  'src/CDT/Triangulation.h'". This file is part of the bundled, unmodified
  'CDT' C++ library; the pragma is upstream and only silences a compiler
  warning in third-party code.

## Bundled code

The package includes the header-only 'CDT' C++ library
(https://github.com/artem-ogre/CDT), licensed under the Mozilla Public
License 2.0 (see inst/LICENSE.CDT), and the robust geometric predicates of
William C. Lenthe (predicates.h), under a 3-clause BSD license. Copyright
holders are acknowledged in Authors@R and inst/COPYRIGHTS.
