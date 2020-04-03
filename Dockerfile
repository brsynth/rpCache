FROM brsynth/rdkit:debian

RUN python -m pip install redis

COPY rpCache.py /home/
