FROM nfcore/base
LABEL authors="Peter Kruczkiewicz" \
      description="Docker image containing all requirements for peterk87/nf-iav-ont pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-iav-ont-1.0dev/bin:$PATH
