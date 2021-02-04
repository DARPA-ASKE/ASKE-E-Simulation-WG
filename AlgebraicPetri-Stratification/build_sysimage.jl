using PackageCompiler;
using Pkg;

create_sysimage(
    [Symbol(x.name) for x in values(Pkg.dependencies()) if x.is_direct_dep];
    sysimage_path="SysImage.so"
)
