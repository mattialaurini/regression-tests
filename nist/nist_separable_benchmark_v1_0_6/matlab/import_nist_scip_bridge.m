function [module, info] = import_nist_scip_bridge(args)
%IMPORT_NIST_SCIP_BRIDGE Import the package-local Python bridge safely.
arguments
    args.PackageRoot (1,1) string = fileparts(fileparts(mfilename("fullpath")))
    args.Reload (1,1) logical = false
end
expectedApiVersion = "1.0.0";
packageRoot = string(args.PackageRoot);
bridgeDirectory = fullfile(packageRoot,"python");
moduleName = "nist_scip_bridge";
if ~isfolder(bridgeDirectory)
    error("NISTSCIP:BridgeDirectoryNotFound", ...
        "Python bridge directory not found: %s",bridgeDirectory);
end
insert(py.sys.path,int32(0),char(bridgeDirectory));
py.importlib.invalidate_caches();
module = py.importlib.import_module(moduleName);
if args.Reload
    module = py.importlib.reload(module);
end
if ~logical(py.builtins.hasattr(module,"bridge_info_json"))
    error("NISTSCIP:BridgeApiMissing", ...
        "The imported module does not expose bridge_info_json().");
end
info = jsondecode(string(module.bridge_info_json()));
actualVersion = string(info.api_version);
actualFile = string(info.module_file);
actualDirectory = string(fileparts(actualFile));
if actualVersion ~= expectedApiVersion
    error("NISTSCIP:BridgeApiMismatch", ...
        "Expected bridge API %s, got %s from %s.", ...
        expectedApiVersion,actualVersion,actualFile);
end
if actualDirectory ~= bridgeDirectory
    error("NISTSCIP:BridgePathMismatch", ...
        ["Python imported the wrong bridge.\nExpected: %s\nImported: %s\n" + ...
         "Restart the out-of-process Python worker and retry."], ...
        bridgeDirectory,actualFile);
end
end
