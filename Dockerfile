# 酶搜索 + DLKcat 流水线镜像（CPU）
FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive \
    PIPELINE_ROOT=/app \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/conda/envs/gmas_pipeline/bin:$PATH

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
        git unzip ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /app/environment.yml
RUN conda config --set channel_priority flexible \
    && conda env create -n gmas_pipeline -f /app/environment.yml \
    && conda clean -afy

# 代码与脚本（不含大体量 results/structures）
COPY src /app/src
COPY scripts /app/scripts
COPY config.yaml /app/config.yaml
COPY config.example.yaml /app/config.example.yaml
COPY docs /app/docs
COPY README.md /app/README.md
COPY requirements.txt /app/requirements.txt
COPY third_party /app/third_party

RUN chmod +x /app/scripts/*.sh \
    && bash /app/scripts/vendor_dlkcat.sh \
    && mkdir -p /app/data /app/results

ENTRYPOINT ["/app/scripts/entrypoint.sh"]
CMD []
