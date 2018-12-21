FROM ubuntu:bionic

LABEL maintainer="yx2@sanger.ac.uk"

USER root

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

RUN mkdir -p /opt
COPY build/build_base.sh /opt
WORKDIR /opt
RUN ./build_base.sh && rm build_base.sh
COPY scripts /opt
RUN chmod a+rx -R /opt

USER ubuntu
ENV PATH "$PATH:/opt"
WORKDIR /home/ubuntu

CMD ["/bin/bash"]