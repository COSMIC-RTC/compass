-- -*- lua -*-
--[[
   COSMIC modulefile

   .modulefiles/cosmic/local.lua
--]]
family("compass")
local module_version = myModuleVersion()

help([[
        io.stdout:write("\tThis module defines the compass environment")
]])

whatis("This module defines the compass environment")

function pushenvVar(var, value)
    if (mode() == "load") then
        pushenv(var, value)
        io.stderr:write("\tUsing " .. var .. "=" .. value .. "\n")
    end
end
local home_local = os.getenv("HOME") .. "/local/"
local compass = os.getenv("HOME") .. "/compass/"
local compass_install = os.getenv("HOME") .. "/local/compass"
local shesha = os.getenv("HOME") .. "/compass/shesha"
pushenv("COMPASS_ROOT", compass)
pushenv("COMPASS_INSTALL_ROOT", compass_install)
pushenv("SHESHA_ROOT", shesha)
prepend_path("PYTHONPATH", shesha)
prepend_path("PYTHONPATH", compass_install .. "/python")
prepend_path("LD_LIBRARY_PATH", compass_install .. "/lib")
prepend_path("PKG_CONFIG_PATH", compass_install .. "/lib/pkgconfig")

if (not isloaded("cuda")) then
    load("cuda")
end

