let (|>) a b = b a
let car, cdr = fst, snd
let odd n = n land 1 = 1
let even n = not (odd n)
module ExtArray = ExtArray.Array
module ExtList = ExtList.List
module ExtString = ExtString.String
module Program = struct
	let name = "Naive Bayes Classifier - Tabulate"
	let version = "1.0"
end
let chop_extra s = ExtString.strip ~chars:"_- " s
let any_alphanumeric s =
	let e = String.length s in
	let rec loop i =
		if i = e then false
		else match s.[i] with
			'a'..'z' | 'A'..'Z' | '0'..'9' -> true
			| _ -> loop (i + 1)
	in loop 0
let guess_prefix filenames =
	let s =
		List.map Misc.basename_without_extension filenames
			|> Misc.longest_common_substring |> chop_extra
	in
	if any_alphanumeric s then s
	else ""
module Options = struct
	let parser = OptParse.OptParser.make ~version:Program.version ()
	let genomes =
		let option = OptParse.StdOpt.str_option ~metavar:"file" () in
		OptParse.OptParser.add parser ~short_name:'g' ~long_name:"genome-list" option;
		(fun () -> match OptParse.Opt.opt option with
			Some x -> x
			| None ->
				OptParse.OptParser.usage parser ();
				exit 1
		)
	let columns =
		let option = OptParse.StdOpt.int_option ~default:200 () in
		OptParse.OptParser.add parser ~short_name:'c' ~long_name:"columns" option;
		(fun () -> OptParse.Opt.get option)
	let template =
		let option = OptParse.StdOpt.str_option ~metavar:"template"
			~default:"$prefix-$n.csv.gz" () in
		OptParse.OptParser.add parser ~short_name:'t' ~long_name:"output-template" option;
		(fun () -> OptParse.Opt.get option)
	let given_prefix =
		let x = ref None in
		OptParse.StdOpt.str_callback ~metavar:"prefix" (fun s -> x := Some s)
			|> OptParse.OptParser.add parser ~short_name:'p'
				~long_name:"output-prefix";
		(fun () -> !x)
	let files = OptParse.OptParser.parse_argv parser
	let prefix = match given_prefix () with
		None -> (
			match guess_prefix files with
				"" -> "table"
				| s -> s
		) | Some x -> x
	let string_of_int n = Printf.sprintf "%08u" n
	let output n = template () |> Misc.substitute ["prefix", prefix; "n", string_of_int n]
end

exception Collide of string * string * string
let match_files_to_genomes genomes files =
	let g = Array.of_list genomes in
	let gl = Array.length g in
	let gi = Array.init gl (fun i -> i) in
	Array.sort (fun ai bi -> compare (String.length g.(bi)) (String.length g.(ai))) gi;
	let gf = Array.make gl None in
	List.iter (fun file ->
		try (
			let i = ExtArray.find (fun i -> ExtString.exists file g.(i)) gi in
			match gf.(i) with
				None -> gf.(i) <- Some file
				| Some other_file -> raise (Collide (g.(i), other_file, file))
		) with Not_found -> () (* file does not match any genomes *)
	) files;
	let r = ref [] in
	for i = gl - 1 downto 0 do
		match gf.(i) with
			None -> () (* no file matches a given genome *)
			| Some file -> r := (g.(i), file) :: !r
	done;
	!r

let columns = Options.columns ()
let genomes, files =
	let genomes = Misc.io_of_file (Options.genomes ()) |> Misc.enum_of_lines
		|> ExtList.of_enum in
	let g, f = match_files_to_genomes genomes Options.files |> List.split in
	Array.of_list g,
	f |> List.map (fun x -> x, x |> open_in |> IO.input_channel) |> Array.of_list

let newfile, newline, cell, finish = 
	let i = ref 0 in
	let open_file () =
		i := !i + 1;
		let filename = Options.output !i in
		Gzip.open_out filename
	in
	let c = ref None in
	let force_out () =
		match !c with
			None ->
				let d = open_file () in
				c := Some d;
				d
			| Some d -> d
	in
	let start_of_line = ref true in
	(* newfile *) (fun () ->
		(match !c with
			Some c -> Gzip.close_out c
			| None -> ());
		i := !i + 1;
		c := None
	),
	(* newline *) (fun () ->
		let c = force_out () in
		Gzip.output_char c '\n';
		start_of_line := true
	),
	(* cell *) (fun s ->
		let c = force_out () in
		if not !start_of_line then Gzip.output_char c ',';
		Gzip.output c s 0 (String.length s);
		start_of_line := false
	),
	(* finish *) (fun () ->
		match !c with
			Some c -> Gzip.close_out c
			| None -> ()
	)

let read_two_fields (filename, c) =
		let line = IO.read_line c in
		match Misc.split2 line with
			Some x -> x
			| None ->
				Printf.eprintf "\
					There is something wrong with the file %s. \
					The offending line is:\n%s\n" filename line;
				exit 1

let output_file () =
	newfile ();
	let yes_we_have_input =
		let x = ref false in
		(fun () ->
			if !x then ()
			else (
				cell "names";
				x := true
			)
		)
	in
	let a = Array.make columns "" in
	let rec loop i =
		if i < columns then (
			try (
				let name, datum = read_two_fields files.(0) in
				yes_we_have_input ();
				cell name;
				a.(i) <- datum;
				loop (i + 1)
			) with IO.No_more_input ->
				if i > 0 then Array.sub a 0 i
				else raise End_of_file
		) else a
	in
	let a = loop 0 in
	let these_columns = Array.length a in
	newline ();
	cell genomes.(0);
	for i = 0 to these_columns - 1 do cell a.(i) done;
	newline ();
	for i = 1 to Array.length genomes - 1 do
		cell genomes.(i);
		for j = 0 to these_columns - 1 do
			let _, datum = read_two_fields files.(i) in
			cell datum
		done;
		newline ()
	done

let () =
	try while true do output_file () done
	with End_of_file -> finish ()
