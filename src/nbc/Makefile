MLTON = mlton
MLYACC = mlyacc
MLLEX = mllex
%.grm.sig %.grm.sml: %.grm
	$(MLYACC) $^
%.lex.sml: %.lex
	$(MLLEX) $^
%: %.mlb
	$(MLTON) $(MLTONFLAGS) -output $@ $^
all: count probabilities-by-read
count: count.mlb
probabilities-by-read: probabilities-by-read.mlb
score: score.mlb
tabulate: tabulate.mlb
clean:
	rm -f count probabilities-by-read
