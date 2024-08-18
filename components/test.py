import gdsfactory as gf
import ubcpdk

@gf.cell
def test():
    c = gf.Component()
    c << ubcpdk.components.dbg()

    return c

if __name__ == "__main__":
    test().show()