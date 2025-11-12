-- -*- lua -*-
--[[
   CUDA modulefile

   .modulefiles/cuda/system.lua
--]]
family("CUDA")
local module_version = myModuleVersion()

help([[
        io.stdout:write("\tThis module defines the CUDA environment")
]])

whatis("This module defines the CUDA environment")

function pushenvVar(var, value)
    if (mode() == "load") then
        pushenv(var, value)
        io.stderr:write("\tUsing " .. var .. "=" .. value .. "\n")
    end
end

local cuda_root = "/usr/local/cuda"
-- local	nvidia_ver="361"

-- pushenv("GENCODE", "arch=compute_61,code=sm_61")
pushenv("CUDA_ROOT", cuda_root)
pushenv("CUDA_INC_PATH", cuda_root .. "/include")
pushenv("CUDA_LIB_PATH", cuda_root .. "/lib")
pushenv("CUDA_LIB_PATH_64", cuda_root .. "/lib64")
-- pushenv("CPLUS_INCLUDE_PATH", cuda_root .. "/include")
prepend_path("PATH", cuda_root .. "/bin")
prepend_path("LD_LIBRARY_PATH", cuda_root .. "/lib64:" .. cuda_root .. "/lib")
-- prepend_path("LD_LIBRARY_PATH", "/usr/lib/nvidia-" .. nvidia_ver)
-- prepend_path("LD_LIBRARY_PATH", "/usr/lib32/nvidia-" .. nvidia_ver)

-- set-alias "startCUDA" 	"sudo tee /proc/acpi/bbswitch <<< ON"
-- set-alias "stopCUDA"	"sudo rmmod nvidia_uvm; sudo rmmod nvidia; sudo tee /proc/acpi/bbswitch <<< OFF"
-- set-alias "nvidia-smi"  "optirun /usr/lib/nvidia-361/bin/nvidia-smi"
