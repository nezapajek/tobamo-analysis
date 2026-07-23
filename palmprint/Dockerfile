FROM gcc:7.2
WORKDIR /usr/src
RUN git clone https://github.com/rcedgar/palmscan.git
WORKDIR /usr/src/palmscan/bin
WORKDIR /usr/src/palmscan/src
RUN make ../bin/palmscan2
WORKDIR /usr/src/palmscan/bin
ENTRYPOINT ["../bin/palmscan2"]