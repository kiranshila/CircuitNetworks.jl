# Indexing operations for an N-port network

import Base.getindex, Base.view

# Parameter access and conversion

# We'll have to be careful here being correct and avoiding copies when indexing down the frequency axis