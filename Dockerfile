FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["bash", "-lc"]
WORKDIR /work

# --- minimal tools to fetch Miniconda ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl bzip2 \
 && rm -rf /var/lib/apt/lists/*

# --- install Miniconda (same effect as continuumio/miniconda3) ---
RUN curl -sSLo /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && rm -f /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# --- Python2 env for Mitsuba Compile ---
COPY environment_mitsuba_compile.yml .
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create -n mitsuba0.5python2 -f environment_mitsuba_compile.yml && conda clean -afy
ENV PATH=/opt/conda/envs/mitsuba0.5python2/bin:$PATH

# --- Python3 env for OHD Tutorial ---
COPY environment_tutorial.yml .
RUN conda env create -n mitsubaohd -f environment_tutorial.yml && conda clean -afy
ENV PATH=/opt/conda/envs/mitsubaohd/bin:$PATH

# --- build dependencies (unchanged from your list) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    scons \
    mercurial \
    libilmbase-dev \
    libpng-dev \
    libjpeg-dev \
    libxerces-c-dev \
    libboost-all-dev \
    libopenexr-dev \
    libglewmx-dev \
    libxxf86vm-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    libeigen3-dev \
    libfftw3-dev \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate mitsuba0.5python2' >> /root/.bashrc
CMD ["/bin/bash","-l"]