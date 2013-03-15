type svalue = Tokens.svalue
type ('a, 'b) token = ('a, 'b) Tokens.token
type pos = int
type lexresult = (svalue, pos) token
type arg = int ref
fun eof _ = Tokens.EOF (~1, ~1)
%%
%full
%header (functor SubstitutionLexFun (structure Tokens: SubstitutionGrm_TOKENS));
%arg (nesting);
%s VARIABLE PARENTHESIZED BRACED;
text = ([^$] | "\\$")+;
dollar = "$";
leftparenthesis = "(";
rightparenthesis = ")";
notparenthesis = [^()]+;
leftbrace = "{";
rightbrace = "}";
notbrace = [^{}]+;
other = [A-Za-z0-9_]+;
%%
<INITIAL>{text} => (Tokens.TEXT (yytext, yypos, yypos + size yytext));
<INITIAL>{dollar} => (YYBEGIN VARIABLE; Tokens.DOLLAR (yypos, yypos + 1));
<INITIAL>{dollar}{leftparenthesis} => (
	YYBEGIN PARENTHESIZED
	; nesting := !nesting + 1
	; Tokens.LEFT_PARENTHESIS (yypos, yypos + 1)
);
<INITIAL>{dollar}{leftbrace} => (
	YYBEGIN BRACED
	; nesting := !nesting + 1
	; Tokens.LEFT_BRACE (yypos, yypos + 1)
);
<VARIABLE>{other} => (YYBEGIN INITIAL; Tokens.TEXT (yytext, yypos, yypos + size yytext));
<PARENTHESIZED>{notparenthesis} => (Tokens.TEXT (yytext, yypos, yypos + size yytext));
<PARENTHESIZED>{leftparenthesis} => (
	nesting := !nesting + 1
	; Tokens.LEFT_PARENTHESIS (yypos, yypos + 1)
);
<PARENTHESIZED>{rightparenthesis} => (
	nesting := !nesting - 1
	; if !nesting = 0 then YYBEGIN INITIAL else ()
	; Tokens.RIGHT_PARENTHESIS (yypos, yypos + 1)
);
<BRACED>{notbrace} => (Tokens.TEXT (yytext, yypos, yypos + size yytext));
<BRACED>{leftbrace} => (
	nesting := !nesting + 1
	; Tokens.LEFT_BRACE (yypos, yypos + 1)
);
<BRACED>{rightbrace} => (
	nesting := !nesting - 1
	; if !nesting = 0 then YYBEGIN INITIAL else ()
	; Tokens.RIGHT_BRACE (yypos, yypos + 1)
);
