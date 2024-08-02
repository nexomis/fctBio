# 1. Installation

First clone the repository to `PATH` and then install on R.

## 1.1 Build documentation

```
apptainer run docker://quay.io/nexomis/r-nexoverse:4.3.3-Bioc_3.18-06.2024 Rscript .ci/build_documentation.r
docker run -u $UID:$GID -v $PWD:$PWD -w $PWD quay.io/nexomis/r-nexoverse:4.4.1-bioc_3.19-07.2024 Rscript .ci/test.r

docker run -u $UID:$GID -v $PWD:$PWD -w $PWD quay.io/nexomis/r-nexoverse:4.4.1-bioc_3.19-07.2024 Rscript .ci/build_documentation.r

```

## 1.2 Build and install the library

```
devtools::install_local(path=PATH,force=T,upgrade = 'never')
```

# Getting started 

See the vignette `Examples`


# check dependencies 

```
cat R/* | grep "::" | sed -re "s/.*[\(\[ ]([a-zA-Z0-9_]+)\:\:.*/\1/g" | sort -u | xargs | sed 's/ /, /g'
```
