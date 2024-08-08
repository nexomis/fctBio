FROM quay.io/nexomis/r-nexoverse:4.4.1-bioc_3.19-07.2024

RUN export DEBIAN_FRONTEND=noninteractive \ 
  && apt-get update \
  && apt-get install -y pandoc

RUN Rscript -e "\
  rspm::enable();\
  install.packages('pandoc', clean = TRUE);\
  install.packages('kableExtra', clean = TRUE)\
"
