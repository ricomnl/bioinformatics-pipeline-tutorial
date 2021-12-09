FROM ubuntu:20.10
RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install matplotlib