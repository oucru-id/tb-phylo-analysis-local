import sys
import importlib

def main():
    packages = [
        "pandas", "networkx", "pyvis", "matplotlib", 
        "seaborn", "numpy", "Bio"
    ]

    print("python_tools:")
    print(f"  python: {sys.version.split()[0]}")
    
    for p in packages:
        try:
            mod = importlib.import_module(p)
            version = getattr(mod, "__version__", "unknown")
            print(f"  {p}: {version}")
        except ImportError:
            print(f"  {p}: not_installed")

if __name__ == "__main__":
    main()