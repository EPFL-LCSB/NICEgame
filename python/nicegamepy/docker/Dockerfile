


FROM python:3.6

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
        less            \
    && rm -rf /var/lib/apt/lists/*

ENV USER nicegamepy
ENV HOME /home/$USER

RUN useradd -ms "/bin/bash" "$USER"
USER $USER
WORKDIR $HOME

USER root

# Copy python package requirements
COPY requirements.txt .

# Install python packages
RUN pip install -r requirements.txt

# Take care of the solvers
COPY ./solvers /solvers
COPY ./utils /utils

RUN chmod u+x /utils/*.sh

# Install CPLEX
RUN /utils/install_cplex.sh

# Remove installers
RUN rm -rf /solvers


RUN mkdir /nicegamepy

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"

#Finalizing installation


USER $USER
RUN mkdir ./work
WORKDIR ./work


# Load your package in development mode on startup
ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /nicegamepy/docker/pytfa && \
pip install --user -e /nicegamepy && \
$0 $*"]

CMD /bin/bash
