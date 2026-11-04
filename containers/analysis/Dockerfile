FROM mambaorg/micromamba:2.0.5
USER root
RUN apt-get update && apt-get install -y --no-install-recommends g++ procps && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge python=3.11.11 pip=25.0.1 ncurses=6.5 && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt && python -m pip check
ARG MODEL_URL=https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl
ARG MODEL_SHA256=290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502
RUN mkdir -p /opt/conda/models && python -c "import urllib.request,hashlib; p='/opt/conda/models/Immune_All_Low.pkl'; urllib.request.urlretrieve('${MODEL_URL}',p); assert hashlib.sha256(open(p,'rb').read()).hexdigest() == '${MODEL_SHA256}'"
ENV PATH=/opt/conda/bin:$PATH MPLBACKEND=Agg MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
USER root
RUN python -c "import scanpy,celltypist,scrublet; from celltypist.models import Model; Model.load('/opt/conda/models/Immune_All_Low.pkl')"
ENTRYPOINT []
CMD ["bash"]
