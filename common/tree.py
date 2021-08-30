"""Collection of tree data structures."""


from typing import Tuple, List


class RBNode(object):
    """Node for left-leaning red-black tree."""

    def __init__(self, k: int, v: str = None) -> None:
        """Create node with key and value."""
        self.key: int = k
        self.val: str = v
        self.red: bool = True
        self.left: RBNode = None
        self.right: RBNode = None

    def compare(self, k: int) -> int:
        if self.key == k:
            return 0
        elif self.key > k:
            return 1
        else:
            return -1


class RBTree(object):
    """Implementation of left-leaning red-black tree data structure."""

    def __init__(self) -> None:
        """Create empty tree."""
        self.root = None

    def inorder(self) -> Tuple[int, str]:
        return self.__inorder2(self.root)

    def __inorder2(self, n: RBNode) -> Tuple[int, str]:
        if n:
            yield from self.__inorder2(n.left)
            yield n.key, n.val
            yield from self.__inorder2(n.right)

    def size(self) -> int:
        return self.__size2(self.root)

    def __size2(self, n: RBNode) -> int:
        return 1 + self.__size2(n.left) + self.__size2(n.right) if n else 0

    def __is_red(self, n: RBNode) -> bool:
        return n.red if n else False

    def search(self, k: int) -> str:
        return self.__search2(k, self.root)

    def __search2(self, k: int, n: RBNode) -> str:
        if n:
            cmp_res: int = n.compare(k)
            if cmp_res > 0:
                return self.__search2(k, n.left)
            elif cmp_res < 0:
                return self.__search2(k, n.right)
            else:
                return n.val
        return None

    def grab_node(self, k: int) -> RBNode:
        return self.__grab_node2(k, self.root)

    def __grab_node2(self, k: int, n: RBNode) -> RBNode:
        if n:
            cmp_res: int = n.compare(k)
            if cmp_res > 0:
                return self.__grab_node2(k, n.left)
            elif cmp_res < 0:
                return self.__grab_node2(k, n.right)
            else:
                return n
        return None

    def __rotate_left(self, n: RBNode) -> RBNode:
        x: RBNode = n.right
        n.right = x.left
        x.left = n
        x.red = n.red
        n.red = True
        return x

    def __rotate_right(self, n: RBNode) -> RBNode:
        x: RBNode = n.left
        n.left = x.right
        x.right = n
        x.red = n.red
        n.red = True
        return x

    def __flip_colors(self, n: RBNode) -> None:
        n.red = not n.red
        n.left.red = not n.left.red
        n.right.red = not n.right.red

    def insert(self, k: int, v: str = None) -> None:
        self.root = self.__insert2(self.root, k, v)
        self.root.red = False

    def __insert2(self, n: RBNode, k: int, v: str = None) -> RBNode:
        if n is None:
            return RBNode(k, v)
        cmp_res: int = n.compare(k)
        if cmp_res > 0:
            n.left = self.__insert2(n.left, k, v)
        elif cmp_res < 0:
            n.right = self.__insert2(n.right, k, v)
        else:
            n.val = v
        if not self.__is_red(n.left) and self.__is_red(n.right):
            n = self.__rotate_left(n)
        if self.__is_red(n.left) and self.__is_red(n.left.left):
            n = self.__rotate_right(n)
        if self.__is_red(n.left) and self.__is_red(n.right):
            self.__flip_colors(n)
        return n


class RBINode(object):
    """Node for left-leaning red-black interval tree."""

    def __init__(self, ks: int, v: str, ke: int = None) -> None:
        """Create node with key and value."""
        self.key: Tuple[int, int] = (ks, ke) if ke else (ks, ks)
        self.max: int = ke if ke else ks
        self.val: str = v
        self.red: bool = True
        self.left: RBINode = None
        self.right: RBINode = None

    def compare(self, k: int) -> int:
        return self.key[0] - k

    def recal_max(self) -> None:
        max_cand: List[int] = [self.key[1]]
        if self.left:
            max_cand.append(self.left.max)
        if self.right:
            max_cand.append(self.right.max)
        self.max = max(max_cand)


class RBITree(object):
    """Implementation of left-leaning red-black interval tree data structure."""

    def __init__(self) -> None:
        """Create empty tree."""
        self.root = None

    def inorder(self) -> Tuple[int, str]:
        return self.__inorder2(self.root)

    def __inorder2(self, n: RBINode) -> Tuple[int, str]:
        if n:
            yield from self.__inorder2(n.left)
            for i in range(n.key[0], n.key[1] + 1):
                yield i, n.val
            yield from self.__inorder2(n.right)

    def size(self) -> int:
        return self.__size2(self.root)

    def __size2(self, n: RBINode) -> int:
        return 1 + self.__size2(n.left) + self.__size2(n.right) if n else 0

    def __is_red(self, n: RBINode) -> bool:
        return n.red if n else False

    def search(self, k: int) -> str:
        return self.__search2(k, self.root)

    def __search2(self, k: int, n: RBINode) -> str:
        if n:
            if k < n.key[0]:
                return self.__search2(k, n.left)
            elif k > n.key[1]:
                if n.left and n.left.max >= k:
                    return self.__search2(k, n.left)
                elif n.right:
                    return self.__search2(k, n.right)
            else:
                return n.val
        return None

    def grab_node(self, k: int) -> RBINode:
        return self.__grab_node2(k, self.root)

    def __grab_node2(self, k: int, n: RBINode) -> RBINode:
        if n:
            if k < n.key[0]:
                return self.__grab_node2(k, n.left)
            elif k > n.key[1]:
                if n.left and n.left.max >= k:
                    return self.__grab_node2(k, n.left)
                elif n.right:
                    return self.__grab_node2(k, n.right)
            else:
                return n
        return None

    def __rotate_left(self, n: RBINode) -> RBINode:
        x: RBINode = n.right
        n.right = x.left
        x.left = n
        x.red = n.red
        n.red = True
        n.recal_max()
        return x

    def __rotate_right(self, n: RBINode) -> RBINode:
        x: RBINode = n.left
        n.left = x.right
        x.right = n
        x.red = n.red
        n.red = True
        n.recal_max()
        return x

    def __flip_colors(self, n: RBINode) -> None:
        n.red = not n.red
        n.left.red = not n.left.red
        n.right.red = not n.right.red

    def remove_overlap(self, ks: int, ke: int = None) -> None:
        if ke is None or ke == ks:
            overlap_node: RBINode = self.grab_node(ks)
            if overlap_node and overlap_node.key[0] != overlap_node.key[1]:
                if overlap_node.key[0] < ks:
                    overlap_node.key = (overlap_node.key[0], ks - 1)
                elif overlap_node.key[1] > ks:
                    overlap_node.key = (ks + 1, overlap_node.key[1])
                elif ks == overlap_node.key[0] and ks == overlap_node.key[1]:
                    overlap_node.val = None

    def insert(self, ks: int, v: str, ke: int = None) -> None:
        self.remove_overlap(ks, ke)
        self.root = self.__insert2(self.root, ks, v, ke)
        self.root.red = False

    def __insert2(self, n: RBINode, ks: int, v: str, ke: int) -> RBINode:
        if n is None:
            return RBINode(ks, v, ke)
        keyCmp: int = n.compare(ks)
        if keyCmp > 0:
            n.left = self.__insert2(n.left, ks, v, ke)
        elif keyCmp < 0:
            n.right = self.__insert2(n.right, ks, v, ke)
        else:
            n.val = v
        if not self.__is_red(n.left) and self.__is_red(n.right):
            n = self.__rotate_left(n)
        if self.__is_red(n.left) and self.__is_red(n.left.left):
            n = self.__rotate_right(n)
        if self.__is_red(n.left) and self.__is_red(n.right):
            self.__flip_colors(n)
        n.recal_max()
        return n
