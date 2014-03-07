PREFIX = "/usr/"
all: c

install: c 
	@echo installing executable files to ${DESTDIR}${PREFIX}/bin
	@mkdir -p ${DESTDIR}${PREFIX}/bin
	@cp -vf src/nbc/probabilities-by-read ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	@cp -vf src/nbc/count ${DESTDIR}${PREFIX}/bin/count-kmers
	@cp -vf src/c/quikr_train ${DESTDIR}${PREFIX}/bin/quikr_train
	@cp -vf src/c/quikr ${DESTDIR}${PREFIX}/bin/quikr
	@cp -vf src/c/multifasta_to_otu ${DESTDIR}${PREFIX}/bin/multifasta_to_otu 
	@cp -vf src/python/generate_kmers ${DESTDIR}${PREFIX}/bin/generate_kmers
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/count-kmers
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/quikr
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/quikr_train
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/multifasta_to_otu
	chmod -v 555 ${DESTDIR}${PREFIX}/bin/generate_kmers
	@cp -vf src/c/quikr.1 ${DESTDIR}${PREFIX}/share/man/man1/quikr.1
	@cp -vf src/c/quikr_train.1 ${DESTDIR}${PREFIX}/share/man/man1/quikr_train.1
	@cp -vf src/c/multifasta_to_otu.1 ${DESTDIR}${PREFIX}/share/man/man1/multifasta_to_otu.1

c:
	@echo "building c"
	@cd src/c; make
python:
	@echo "configuring python"
	@cd src/python; python setup.py build

install_python:
	@cd src/python; python setup.py install
	@cp -vf src/nbc/probabilities-by-read ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	@cp -vf src/nbc/count ${DESTDIR}${PREFIX}/bin/count-kmers
	@cp -vf src/python/quikr ${DESTDIR}${PREFIX}/bin/quikr.py
	@cp -vf src/python/quikr_train ${DESTDIR}${PREFIX}/bin/quikr_train.py
	@cp -vf src/python/multifasta_to_otu ${DESTDIR}${PREFIX}/bin/multifasta_to_otu.py
	@cp -vf src/python/generate_kmers ${DESTDIR}${PREFIX}/bin/generate_kmers
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/count-kmers
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/quikr.py
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/quikr_train.py
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/multifasta_to_otu.py
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/generate_kmers

clean:
	@echo "cleaning up"
	@cd src/python; rm build -Rvf
	@cd src/nbc; make clean
	@cd src/c; make clean
