const table="accession2taxid"

function create!(source::String, db::SQLite.DB ; chunk::Int=400, overwrite::Bool=false, header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    if overwrite
        DBInterface.execute(db,"DROP TABLE IF EXISTS $table")
    end
    DBInterface.execute(db,"CREATE TABLE $table(accession TEXT PRIMARY KEY, taxid INTEGER)")
    insert!(db, source; chunk=chunk, header=header, accession_col=accession_col, taxid_col=taxid_col)
    return db
end

function Base.insert!(db::SQLite.DB, source::String; chunk::Int=chunk, header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    f = open(source, "r")
    if header
        readline(f) #skip header
    end
    _insert_rows!(db, f, chunk, accession_col, taxid_col)
    close(f)
end

function _insert_rows!(db::SQLite.DB, f::IOStream, chunk::Int, accession_col::Int, taxid_col::Int)
    while !eof(f)
        chunks = [readline(f) for _ in 1:chunk]
        cols = map( x-> split(x, "\t"), chunks)
        accession2taxid = map(x -> (getindex(x,accession_col),getindex(x,taxid_cols), cols)
        _accession2taxid = map(x -> "("*first(x)*","*last(x)*")", accession2taxid))
        rows = join(_accession2taxid)
        DBInterface.execute(db, "INSERT INTO $table(accession, taxid) VALUES $rows")
    end
end

function Base.get(db::SQLite.DB, accession::String, default::Any)
    result = DBInterface.execute(db, "SELECT * FROM $table WHERE accession = \"$accession\"")
    id = [row.taxid for row in result]
    @assert length(id) < 2
    length(id) == 0 && return default
    return id[1]
end