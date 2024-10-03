library(future)
library(RSQLite)

# Plan for parallel processing (modify as needed)
plan(multisession)

# Connect to the database
db <- dbConnect(SQLite(), dbname = "chrom3_1.db")

# Create the table
dbExecute(db, "CREATE TABLE IF NOT EXISTS chrom_data (identifier TEXT PRIMARY KEY, info TEXT)")

# Open the file connection in binary mode
con <- file("data/input/annotations/Chrom3_1.trn", "rb")

# Use a buffer for reading
buffer <- raw(1000)  # Adjust buffer size if needed

while(TRUE) {
  # Read a chunk of bytes into the buffer
  n <- readBin(con, raw(), n = length(buffer))
  if (n[1] == 0) break  # Check for end of file
  
  # Convert bytes to characters, ensuring raw vector input
  lines <- rawToChar(buffer[1:n[1]], multiple = TRUE)  # Specify multiple = TRUE
  lines <- strsplit(lines, "\n", fixed = TRUE)[[1]]
  
  # Process and insert lines in parallel
  future({
    dbBeginTransaction(db)
    for (line in lines) {
      if (nchar(line) > 0) {
        id <- substr(line, 1, 10)
        info <- substr(line, 12)
        dbExecute(db, "INSERT INTO chrom_data VALUES (?, ?)", params = list(id, info))
      }
    }
    dbCommit(db)
  }) %plan% multisession
}

# Close connections
dbDisconnect(db)
close(con)