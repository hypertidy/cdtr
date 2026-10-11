## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.


## Bundled code

The package includes the header-only 'CDT' C++ library
(https://github.com/artem-ogre/CDT), licensed under the Mozilla Public
License 2.0 (see inst/LICENSE.CDT), and the robust geometric predicates of
William C. Lenthe (predicates.h), under a 3-clause BSD license. Copyright
holders are acknowledged in Authors@R and inst/COPYRIGHTS.

The bundled 'src/CDT/Triangulation.h' has been patched to remove a pragma
that suppressed '-Wunused-parameter' diagnostics (the unused parameter names
are now commented out instead); the change is documented in inst/COPYRIGHTS.
