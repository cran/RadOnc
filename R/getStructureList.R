setGeneric("getStructureList",
	function (x, ...) {
		standardGeneric("getStructureList")	
	}
)

setMethod("getStructureList", "DVH",
	function (x, structures) {
		as(x, "DVH.list")			
	}
)

setMethod("getStructureList", "DVH.list",
	function (x, structures) {
		return(x[structures])			
	}
)

setMethod("getStructureList", "list",
	function (x, structures) {
		as(lapply(x, getStructureList, structures=structures), "DVH.list")	
	}
)
