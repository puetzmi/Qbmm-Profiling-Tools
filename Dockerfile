FROM intel/oneapi-basekit:2022.2-devel-ubuntu20.04

COPY . /qmom-performance-profiling

RUN apt update
RUN apt install libeigen3-dev -y

RUN pip install quadmompy
