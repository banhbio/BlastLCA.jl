const table="accession2taxid"

function create!(source::String, db::SQLite.DB ; overwrite::Bool=false, header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    if overwrite
        DBInterface.execute(db,"DROP TABLE IF EXISTS $table")
    end
    DBInterface.execute(db,"CREATE TABLE $table(accession TEXT PRIMARY KEY, taxid INTEGER)")
    insert!(db, source; header=header, accession_col=accession_col, taxid_col=taxid_col)
    return db
end

function Base.insert!(db::SQLite.DB, source::String; header::Bool=true, accession_col::Int=2, taxid_col::Int=3)
    f = open(source, "r")
    if header
        readline(f) #skip header
    end
    _insert_rows!(db, f, accession_col, taxid_col)
    close(f)
end

function _insert_rows!(db::SQLite.DB, f::IOStream, accession_col::Int, taxid_col::Int)
    for line in eachline(f)
        cols = split(line, "\t")
        accession = cols[accession_col]
        taxid = cols[taxid_col]
        stmt = DBInterface.prepare(db, "INSERT INTO $table(accession, taxid) VALUES (?, ?)")
        DBInterface.execute(stmt, [accession, taxid])
    end
end

function Base.get(db::SQLite.DB, accession::String, default::Any)
    result = DBInterface.execute(db, "SELECT * FROM $table WHERE accession = \"$accession\"")
    id = [row.taxid for row in result]
    @assert length(id) < 2
    length(id) == 0 && return default
    return id[1]
end