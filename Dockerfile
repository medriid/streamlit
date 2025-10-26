FROM continuumio/miniconda3
WORKDIR /app
COPY environment.yml ./
RUN conda env create -f environment.yml
SHELL ["conda", "run", "-n", "molevis", "/bin/bash", "-lc"]
# Copy app
COPY . /app
# Build tailwind css if tailwind cli is available in container (optional)
# Expose default streamlit port
EXPOSE 8501
ENV PATH /opt/conda/envs/molevis/bin:$PATH
CMD ["conda", "run", "-n", "molevis", "streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]