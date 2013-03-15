all: nbc python

install: 
	@echo installing executable files to ${DESTDIR}${PREFIX}/bin
	@mkdir -p ${DESTDIR}${PREFIX}/bin
	@cp -vf src/nbc/probabilities-by-read ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	@cp -vf src/nbc/count ${DESTDIR}${PREFIX}/bin/count-kmers
	@cp -vf src/python/quikr ${DESTDIR}${PREFIX}/bin/quikr
	@cp -vf src/python/quikr_train ${DESTDIR}${PREFIX}/bin/quikr_train
	@cp -vf src/python/multifasta_to_otu ${DESTDIR}${PREFIX}/bin/multifasta_to_otu
	@cp -vf src/python/generate_kmers ${DESTDIR}${PREFIX}/bin/generate_kmers
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/probabilities-by-read
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/count-kmers
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/quikr
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/quikr_train
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/multifasta_to_otu
	chmod -v 755 ${DESTDIR}${PREFIX}/bin/generate_kmers
	@cd src/python; python setup.py install

nbc:
	@echo "building nbc"
	@cd src/nbc; make

python:
	@echo "configuring python"
	@cd src/python; python setup.py build

clean:
	@echo "cleaning up"
	@cd src/python; rm build -Rvf
	@cd src/nbc; make clean
