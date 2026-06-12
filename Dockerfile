# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: EUPL-1.2
#
# Full reproducibility toolchain for clothoid-halley-coq.
#
# Captures every tool used to produce the committed artefacts:
#   - Coq 8.20.1 + Coquelicot 3.4.3   (verifies coq/*.v)
#   - Python 3 + numpy/scipy/sympy/matplotlib
#                                     (runs python/ + bar charts)
#   - .NET 8 SDK                      (builds + tests csharp/)
#   - OpenJDK 21 (Temurin) + Maven 3  (builds + tests java/)
#   - Node.js 22 + npm                (builds + tests typescript/)
#   - TeX Live (full)                 (rebuilds the paper PDF)
#
# Build:
#     podman build -t clothoid-halley:latest .
#
# Use:
#     podman run --rm -it -v "$PWD:/workspace" -w /workspace \
#         clothoid-halley:latest bash
#
# Inside the container, the end-to-end reproduction is:
#     make -C coq
#     python3 python/build_golden_vectors.py
#     dotnet test csharp/Clothoid.Halley.Tests
#     mvn -f java/pom.xml test
#     ( cd typescript && npm install && npm test )
#     python3 python/run_all_benches.py
#     python3 docs/mathematics/generate_benchmark_graphs.py
#     ( cd docs/mathematics && pdflatex Clothoid_L_Halley_Solver.tex && \
#       bibtex Clothoid_L_Halley_Solver && \
#       pdflatex Clothoid_L_Halley_Solver.tex && \
#       pdflatex Clothoid_L_Halley_Solver.tex )

ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}

# OCI image labels
LABEL org.opencontainers.image.title="clothoid-halley-coq toolchain" \
      org.opencontainers.image.description="Coq 8.20.1, Python, .NET 8, JDK 21, Node 22, TeX Live for clothoid-halley-coq reproduction" \
      org.opencontainers.image.source="https://github.com/grootstebozewolf/clothoid-halley-coq" \
      org.opencontainers.image.licenses="EUPL-1.2" \
      org.opencontainers.image.authors="Jeroen Bloemscheer <jeroen.bloemscheer@merkator.com>"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Europe/Amsterdam \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ------------------------------------------------------------------
# 1. Base system + build tools + TeX Live (single apt invocation
#    minimises layer count and final image size)
# ------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget gnupg2 git make sudo \
        build-essential pkg-config m4 libgmp-dev \
        unzip xz-utils tzdata locales \
        python3 python3-pip python3-venv \
        opam \
        texlive-latex-recommended \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-fonts-extra \
        texlive-science \
        texlive-bibtex-extra \
        biber \
    && rm -rf /var/lib/apt/lists/*

# Python scientific stack (matches what the scripts import)
RUN pip3 install --no-cache-dir --break-system-packages \
        numpy scipy sympy matplotlib

# ------------------------------------------------------------------
# 2. Coq 8.20.1 + Coquelicot 3.4.3 (opam-built; longest single step)
# ------------------------------------------------------------------
RUN opam init --bare --disable-sandboxing --auto-setup --yes && \
    opam switch create coq8.20 ocaml-base-compiler.4.14.2 --yes && \
    opam repo add coq-released https://coq.inria.fr/opam/released --switch=coq8.20 && \
    opam install -y --switch=coq8.20 coq.8.20.1 coq-coquelicot.3.4.3 && \
    opam clean -a -c -s --logs
# Make opam env active in every shell
RUN echo 'eval $(opam env --switch=coq8.20 --set-switch)' > /etc/profile.d/opam.sh && \
    chmod +x /etc/profile.d/opam.sh
ENV PATH=/root/.opam/coq8.20/bin:${PATH} \
    OCAML_TOPLEVEL_PATH=/root/.opam/coq8.20/lib/toplevel \
    CAML_LD_LIBRARY_PATH=/root/.opam/coq8.20/lib/stublibs:/root/.opam/coq8.20/lib/ocaml/stublibs

# ------------------------------------------------------------------
# 3. .NET 8 SDK
# ------------------------------------------------------------------
RUN curl -fsSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh && \
    chmod +x /tmp/dotnet-install.sh && \
    /tmp/dotnet-install.sh --channel 8.0 --install-dir /usr/share/dotnet && \
    ln -s /usr/share/dotnet/dotnet /usr/local/bin/dotnet && \
    rm /tmp/dotnet-install.sh
ENV DOTNET_ROOT=/usr/share/dotnet \
    DOTNET_NOLOGO=1 \
    DOTNET_CLI_TELEMETRY_OPTOUT=1

# ------------------------------------------------------------------
# 4. OpenJDK 21 (Eclipse Temurin) + Maven 3
# ------------------------------------------------------------------
RUN mkdir -p /etc/apt/keyrings && \
    wget -qO - https://packages.adoptium.net/artifactory/api/gpg/key/public \
        | gpg --dearmor -o /etc/apt/keyrings/adoptium.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/adoptium.gpg] https://packages.adoptium.net/artifactory/deb $(. /etc/os-release && echo $VERSION_CODENAME) main" \
        > /etc/apt/sources.list.d/adoptium.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends temurin-21-jdk maven && \
    rm -rf /var/lib/apt/lists/*
# JAVA_HOME picks the architecture-appropriate path at runtime
RUN JDK_PATH=$(dirname $(dirname $(readlink -f $(which javac)))) && \
    echo "export JAVA_HOME=${JDK_PATH}" > /etc/profile.d/java.sh && \
    chmod +x /etc/profile.d/java.sh
ENV JAVA_HOME=/usr/lib/jvm/temurin-21-jdk-amd64

# ------------------------------------------------------------------
# 5. Node.js 22 + npm (NodeSource)
# ------------------------------------------------------------------
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash - && \
    apt-get install -y --no-install-recommends nodejs && \
    rm -rf /var/lib/apt/lists/* && \
    npm config set update-notifier false

# ------------------------------------------------------------------
# 6. Verification: every toolchain advertises its version at build time
#    (fails the image build if any install step silently broke)
# ------------------------------------------------------------------
RUN bash -lc 'set -e && \
    echo "--- toolchain self-test ---" && \
    coqc --version && \
    python3 --version && \
    dotnet --version && \
    java -version && \
    mvn -v | head -1 && \
    node --version && \
    npm --version && \
    pdflatex --version | head -1 && \
    bibtex --version | head -1 && \
    echo "--- ok ---"'

WORKDIR /workspace
SHELL ["/bin/bash", "-lc"]
CMD ["bash"]
