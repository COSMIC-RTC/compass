sed -i -e 's:123ACHANGER123:<script src="install_script.js"></script>:' public/md_doc_html_install_script.html
sed -i -e 's@compass/commits/master">@compass/commits/master"><img alt="Master status" src="https://gitlab.obspm.fr/compass/compass/badges/master/pipeline.svg"/>@' public/index.html
sed -i -e 's@compass/commits/develop">@compass/commits/develop"><img alt="Develop status" src="https://gitlab.obspm.fr/compass/compass/badges/develop/pipeline.svg"/>@' public/index.html
sed -i -e 's@coverage/index.html">@coverage/index.html"><img alt="coverage report" src="https://gitlab.obspm.fr/compass/compass/badges/develop/coverage.svg"/>@' public/index.html
# sed -i -e 's@<object type="image/svg+xml" data="https://gitlab.obspm.fr/compass/compass/badges/develop/coverage.svg" style="pointer-events: none;">coverage report</object>@<img alt="Develop coverage" src="https://gitlab.obspm.fr/compass/compass/badges/develop/coverage.svg"/>@' public/index.html
