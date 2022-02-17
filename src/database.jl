const table="accession2taxid"

function create!(db::SQLite.DB, source::String; overwrite::Bool=false, header::Bool=true, delim::Union{Nothing, Char, String}=nothing , accession_col::Int=2, taxid_col::Int=3)
    overwrite ? DBInterface.execute(db,"DROP TABLE IF EXISTS $table") : nothing
    col_width = max(accession_col, taxid_col)
    temp_colname= [ "Column_$i" for i in 1:col_width ]
    temp_colname[accession_col] = "accession"
    temp_colname[taxid_col] = "taxid"
    CSV.File(source; header=temp_colname, select=[accession_col,taxid_col], types=Dict(accession_col => String, taxid_col => Int), delim=delim) |> SQLite.load!(db, table)
    SQLite.createindex!(db, table, "accessionindex", "accession")
    return db
end

function Base.get(db::SQLite.DB, accession::String, default::Any)
    result = DBInterface.execute(db, "SELECT * FROM $table WHERE accession = \"$accession\"")
    id = [row.taxid for row in result]
    @assert length(id) < 2
    length(id) == 0 && return default
    return id[1]
end