const table="accession2taxid"

function create!(db::SQLite.DB, source::String; overwrite::Bool=false, header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    overwrite ? DBInterface.execute(db,"DROP TABLE IF EXISTS $table") : nothing
    DBInterface.execute(db,"CREATE TABLE $table(accession TEXT PRIMARY KEY, taxid INTEGER)")
    insert!(db, source; header=header, accession_col=accession_col, taxid_col=taxid_col)
    return db
end

function Base.insert!(db::SQLite.DB, source::String; header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    f = open(source, "r")
    header ? readline(f) : nothing #skip header
    _insert_rows!(db, f, accession_col, taxid_col)
    close(f)
end

function _insert_rows!(db::SQLite.DB, f::IOStream, accession_col::Int, taxid_col::Int)
    lines = [line for line in readlines(f)]
    cols = map( x-> split(x, "\t"), chunks)
    accessions = map(x -> getindex(x,accession_col), cols)
    taxids = map(x -> getindex(x,taxid_col), cols) 
    stmt = DBInterface.prepare(db, "INSERT INTO $table(accession, taxid) VALUES(:accession,:taxid)")
    DBInterface.executemany(stmt, (accession=accessions, taxid=taxids))
end

function Base.get(db::SQLite.DB, accession::String, default::Any)
    result = DBInterface.execute(db, "SELECT * FROM $table WHERE accession = \"$accession\"")
    id = [row.taxid for row in result]
    @assert length(id) < 2
    length(id) == 0 && return default
    return id[1]
end