module ExtArray = ExtArray.Array
let (|>) a b = b a

module Options = struct
	let parser = OptParse.OptParser.make ~version:"1.0" ()
	let per_word =
		let option = OptParse.StdOpt.str_option ~metavar:"file" () in
		OptParse.OptParser.add parser ~short_name:'w' ~long_name:"per-word"
			~help:"file to store per-word counts in" option;
		(fun () -> match OptParse.Opt.opt option with
			Some x -> x
			| None ->
				OptParse.OptParser.usage parser ();
				exit 1
		)
	let total =
		let option = OptParse.StdOpt.str_option ~metavar:"file" () in
		OptParse.OptParser.add parser ~short_name:'t' ~long_name:"total"
			~help:"file to store total count in" option;
		(fun () -> match OptParse.Opt.opt option with
			Some x -> x
			| None ->
				OptParse.OptParser.usage parser ();
				exit 1
		)
	let order =
		let option = OptParse.StdOpt.int_option ~default:15 () in
		OptParse.OptParser.add parser ~short_name:'r' ~long_name:"order" option;
		(fun () -> OptParse.Opt.get option)
	let files = OptParse.OptParser.parse_argv parser
	let total = total ()
	let per_word = per_word ()
	let order = order ()
end

let load_file order judy total name =
	Misc.io_of_gzip name |> Fasta.enum_words order
	|> Enum.fold (fun word total ->
		Judy.bump judy word;
		Judy.bump judy (Gene.reverse word);
		total + 2
	) total
let load_files order names =
	let judy = Judy.create () in
	List.fold_left (load_file order judy) 0 names, judy
let gzip_output_string c s = Gzip.output c s 0 (String.length s)
let () =
	let total_words, judy = load_files Options.order Options.files in
	(
		let c = open_out Options.total in
		output_string c (string_of_int total_words ^ "\n");
		close_out c
	); (
		let c = Gzip.open_out Options.per_word in
		Judy.iter (fun word count ->
			gzip_output_string c
				(String.concat "" [string_of_int count; " "; word; "\n"])
		) judy;
		Gzip.close_out c
	)
