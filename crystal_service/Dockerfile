FROM continuumio/miniconda3
WORKDIR /app
COPY requirements.txt ./
RUN conda create -n crystal -y python=3.11 && \
    /opt/conda/bin/conda run -n crystal pip install --upgrade pip && \
    /opt/conda/bin/conda run -n crystal pip install -r requirements.txt
COPY . /app
ENV PATH /opt/conda/envs/crystal/bin:$PATH
EXPOSE 8000
CMD ["/opt/conda/bin/conda", "run", "-n", "crystal", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
