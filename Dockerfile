FROM brsynth/rdkit:debian

RUN python -m pip install redis

WORKDIR /home

COPY rpCache.py .

RUN ldconfig
