import gdsfactory as gf
import ubcpdk
gf.clear_cache()

@gf.cell
def test():
    c = gf.Component()
    # c << gf.components.straight(width=0.1)
    c << ubcpdk.components.ring_double()

    return c

if __name__ == "__main__":
    test().show()