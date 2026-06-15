"""
SHA-256 Visualization — manim-slides edition  (complete, polished, LaTeX-migrated)
===================================================================================

Render:
    manim-slides render sha256_presentation.py SHA256Presentation --quality h

Present:
    manim-slides present SHA256Presentation

Export to HTML:
    manim-slides convert SHA256Presentation presentation.html

Render individual scenes for quick testing:
    manim -pql sha256_presentation.py TitleSlide
    manim -pql sha256_presentation.py MessageToBits
    manim -pql sha256_presentation.py Preprocessing
    manim -pql sha256_presentation.py IVAndConstants
    manim -pql sha256_presentation.py MerkleDamgard
    manim -pql sha256_presentation.py MessageSchedule
    manim -pql sha256_presentation.py BitwiseOperations
    manim -pql sha256_presentation.py CompressionRound
    manim -pql sha256_presentation.py DaviesMeyer
    manim -pql sha256_presentation.py AvalancheEffect
    manim -pql sha256_presentation.py FinalSummary

Slide navigation:
    SPACE / →   next slide
    ←           previous slide
    R           replay current slide
    Q           quit
"""

from manim import *
from manim_slides import Slide

# ============================================================================
# 0.  FIPS 180-4 SHA-256 — pure-Python reference implementation
# ============================================================================

MASK32 = 0xFFFFFFFF

K_CONSTANTS = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

H0_INIT = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32


def Sigma0(x):  return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x):  return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sigma0(x):  return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x):  return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(x, y, z):  return ((x & y) ^ (~x & z)) & MASK32
def Maj(x, y, z): return (x & y) ^ (x & z) ^ (y & z)


def pad_message(msg_bytes):
    ell = len(msg_bytes) * 8
    out = msg_bytes + b'\x80'
    while (len(out) * 8) % 512 != 448:
        out += b'\x00'
    out += ell.to_bytes(8, 'big')
    return out


def to_blocks(padded):
    return [padded[i:i + 64] for i in range(0, len(padded), 64)]


def message_schedule(block):
    W = [int.from_bytes(block[i * 4:i * 4 + 4], 'big') for i in range(16)]
    for t in range(16, 64):
        W.append((sigma1(W[t - 2]) + W[t - 7] + sigma0(W[t - 15]) + W[t - 16]) & MASK32)
    return W


def compress_block(block, H_in):
    W = message_schedule(block)
    a, b, c, d, e, f, g, h = H_in
    rounds = []
    for t in range(64):
        T1 = (h + Sigma1(e) + Ch(e, f, g) + K_CONSTANTS[t] + W[t]) & MASK32
        T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
        rounds.append(dict(
            t=t, W=W[t], K=K_CONSTANTS[t],
            before=(a, b, c, d, e, f, g, h),
            T1=T1, T2=T2,
            after=((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g),
        ))
        a, b, c, d, e, f, g, h = (T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g
    H_out = [(x + y) & MASK32 for x, y in zip([a, b, c, d, e, f, g, h], H_in)]
    return H_out, rounds, W


def sha256_full(msg_bytes):
    padded  = pad_message(msg_bytes)
    blocks  = to_blocks(padded)
    H       = H0_INIT[:]
    all_rounds, all_W, all_H_in = [], [], []
    for block in blocks:
        all_H_in.append(H[:])
        H, rounds, W = compress_block(block, H)
        all_rounds.append(rounds)
        all_W.append(W)
    digest = b''.join(x.to_bytes(4, 'big') for x in H)
    return digest, H, padded, blocks, all_rounds, all_W, all_H_in


# ============================================================================
# 1.  COLOUR PALETTE
# ============================================================================

BG       = "#0F0F1A"
AMBER    = "#FFD166"
TEAL     = "#06D6A0"
CORAL    = "#EF476F"
SKY      = "#58C4DD"
LAVENDER = "#B28DFF"
DKGRAY   = "#1E1E2E"
MDGRAY   = "#2A2A3C"
BIT_ON   = SKY
BIT_OFF  = "#2A2A3E"
HIGHLIGHT = AMBER
DIFF_CLR  = CORAL

T1_COLOR  = "#FF8C61"
T2_COLOR  = "#5BC0BE"
NEW_COLOR = AMBER
REG_NAMES = list("abcdefgh")
REG_W, REG_H, REG_BUFF, REG_FONT = 1.58, 1.0, 0.08, 20


# ============================================================================
# 2.  VISUAL HELPERS
# ============================================================================

def byte_to_bits(b):  return format(b, '08b')
def bytes_to_bits(b): return ''.join(byte_to_bits(x) for x in b)
def int_to_bits(v, w):return format(v, f'0{w}b')


def bit_sq(bit_char, size=0.32, on=BIT_ON, off=BIT_OFF):
    col = on if bit_char == '1' else off
    return Square(side_length=size, fill_color=col, fill_opacity=1,
                  stroke_color=WHITE, stroke_width=0.4)


def bit_row(bits, size=0.32, buff=0.04, on=BIT_ON, off=BIT_OFF):
    g = VGroup(*[bit_sq(b, size, on, off) for b in bits])
    g.arrange(RIGHT, buff=buff)
    return g


def brow(bits, sq=0.30, on=SKY, off="#2A2A3E", stroke_w=0.4):
    cells = VGroup(*[
        Square(side_length=sq, fill_color=(on if b == '1' else off),
               fill_opacity=1, stroke_color="#505068", stroke_width=stroke_w)
        for b in bits
    ])
    cells.arrange(RIGHT, buff=0.03)
    return cells


def hex_box(value, width=1.6, height=0.7, font_size=26, color=WHITE):
    box   = RoundedRectangle(corner_radius=0.08, width=width, height=height,
                              stroke_color=color, fill_color=DKGRAY, fill_opacity=1)
    label = Text(f"{value:08x}", font="monospace", font_size=font_size, color=color)
    label.move_to(box.get_center())
    return VGroup(box, label)


def labeled_hex(name, value, width=1.7, height=1.0, nsize=26, vsize=24, nc=AMBER, vc=WHITE):
    box   = RoundedRectangle(corner_radius=0.08, width=width, height=height,
                              stroke_color=nc, fill_color=DKGRAY, fill_opacity=1)
    nt    = Text(name, font_size=nsize, color=nc)
    vt    = Text(f"{value:08x}", font="monospace", font_size=vsize, color=vc)
    nt.move_to(box.get_center() + UP * 0.20)
    vt.move_to(box.get_center() + DOWN * 0.18)
    return VGroup(box, nt, vt)


def make_register(name, value, color=WHITE):
    box   = RoundedRectangle(corner_radius=0.08, width=REG_W, height=REG_H,
                              stroke_color=color, stroke_width=2,
                              fill_color=DKGRAY, fill_opacity=1)
    nt    = Text(name, font_size=18, color=GRAY_B)
    vt    = Text(f"{value:08x}", font="monospace", font_size=REG_FONT)
    nt.move_to(box.get_top()    + DOWN * 0.22)
    vt.move_to(box.get_center() + DOWN * 0.10)
    g = VGroup(box, nt, vt)
    g.box, g.name_t, g.val_t = box, nt, vt
    return g


def make_reg_row(values, y=0.0):
    row = VGroup(*[make_register(n, v) for n, v in zip(REG_NAMES, values)])
    row.arrange(RIGHT, buff=REG_BUFF)
    row.move_to(UP * y)
    return row


def section_title(text, color=WHITE, size=42):
    t = Text(text, font_size=size, color=color)
    t.to_edge(UP, buff=0.40)
    return t


def subtitle_text(text, color=GRAY_B, size=25):
    return Text(text, font_size=size, color=color)


def card_bg(width=11.5, height=1.15):
    return RoundedRectangle(corner_radius=0.14, width=width, height=height,
                             fill_color=MDGRAY, fill_opacity=1,
                             stroke_color="#404060", stroke_width=1.5)


# ============================================================================
# SCENE 0 — Title / Overview
# ============================================================================

class TitleSlide(Slide):
    def construct(self):
        title = Text("SHA-256", font_size=96, color=AMBER, weight=BOLD)
        sub   = Text("Secure Hash Algorithm · 256 bits", font_size=34, color=SKY)
        desc  = Text(
            "A deterministic, one-way function mapping any message\n"
            "to a fixed 256-bit fingerprint",
            font_size=24, color=GRAY_B, line_spacing=1.4, slant=ITALIC)
        VGroup(title, sub, desc).arrange(DOWN, buff=0.4).move_to(ORIGIN + UP * 0.5)

        self.play(Write(title, run_time=1.4))
        self.play(FadeIn(sub, shift=UP * 0.2))
        self.play(FadeIn(desc, shift=UP * 0.15))
        self.next_slide()

        digest = sha256_full(b"abc")[0]
        demo_msg    = Text('"abc"',                 font="monospace", font_size=20, color=TEAL)
        arrow       = Arrow(LEFT, RIGHT, buff=0.1, color=AMBER, stroke_width=3)
        demo_digest = Text(digest.hex(), font="monospace", font_size=12, color=LAVENDER)
        demo_row    = VGroup(demo_msg, arrow, demo_digest).arrange(RIGHT, buff=0.35)
        demo_row.to_edge(DOWN, buff=1.1)

        sha_label = Text("SHA-256", font_size=18, color=AMBER)
        sha_label.next_to(arrow, UP, buff=0.08)

        self.play(Write(demo_msg))
        self.play(GrowArrow(arrow), Write(sha_label))
        self.play(Write(demo_digest))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])
        map_title = Text("What we'll cover", font_size=38, color=WHITE)
        map_title.to_edge(UP, buff=0.5)
        self.play(Write(map_title))

        steps = [
            ("1 · Encoding",         "Message bytes → bits",                 SKY),
            ("2 · Padding",          "Stretch to a multiple of 512 bits",    TEAL),
            ("3 · IV & Constants",   "Where H₀…H₇ and K₀…K₆₃ come from",     AMBER),
            ("4 · Merkle–Damgård",   "How multiple blocks chain together",   LAVENDER),
            ("5 · Message Schedule", "W₀…W₁₅ expanded to W₀…W₆₃",            SKY),
            ("6 · Bitwise Ops",      "Ch, Maj, Σ₀, Σ₁, σ₀, σ₁ explained",    CORAL),
            ("7 · 64-round Loop",    "The compression function step-by-step",AMBER),
            ("8 · Davies–Meyer",     "Feedforward: why it's one-way",        TEAL),
            ("9 · Avalanche",        "1 bit flip → ~128 bits change",        CORAL),
        ]
        
        titles = VGroup()
        descs = VGroup()
        
        for num, desc_str, col in steps:
            titles.add(Text(num, font_size=22, color=col))
            descs.add(Text(desc_str, font_size=20, color=GRAY_B))

        titles.arrange(DOWN, buff=0.28, aligned_edge=LEFT)
        max_title_right = max(t.get_right()[0] for t in titles)
        buffer = 0.8
        
        for t, d in zip(titles, descs):
            d.set_x(max_title_right + buffer, direction=LEFT)
            d.set_y(t.get_y())

        rows = VGroup(*[VGroup(t, d) for t, d in zip(titles, descs)])
        rows.next_to(map_title, DOWN, buff=0.55)
        rows.shift(LEFT * 0.5)

        self.play(LaggedStart(
            *[FadeIn(r, shift=RIGHT * 0.25) for r in rows],
            lag_ratio=0.12
        ))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 1 — Message → bits
# ============================================================================

class MessageToBits(Slide):
    def construct(self):
        message   = "abc"
        msg_bytes = message.encode("ascii")

        title = section_title("Step 1 · Encoding: Message → Bits", color=SKY)
        self.play(Write(title))

        msg_text = Text(f'"{message}"', font="monospace", font_size=72, color=TEAL)
        msg_text.move_to(UP * 0.8)
        self.play(FadeIn(msg_text, shift=UP * 0.3))
        self.next_slide()

        char_group = VGroup()
        for ch in message:
            c_box = VGroup(
                Text(f"'{ch}'", font="monospace", font_size=40, color=TEAL),
                Text(f"ASCII  {ord(ch)}  =  0x{ord(ch):02X}",
                     font_size=21, color=GRAY_B),
            ).arrange(DOWN, buff=0.15)
            char_group.add(c_box)
        char_group.arrange(RIGHT, buff=1.5)
        char_group.next_to(msg_text, DOWN, buff=0.9)

        self.play(FadeOut(msg_text))
        self.play(LaggedStart(*[FadeIn(c, shift=DOWN * 0.3) for c in char_group],
                               lag_ratio=0.2))
        self.next_slide()

        bit_rows = VGroup()
        for ch, c_box in zip(message, char_group):
            bits = byte_to_bits(ord(ch))
            row  = bit_row(bits, size=0.44)
            row.next_to(c_box, DOWN, buff=0.55)
            bit_rows.add(row)

        self.play(*[FadeIn(r, shift=DOWN * 0.2) for r in bit_rows])
        bit_labels = VGroup(*[
            Text(byte_to_bits(ord(ch)), font="monospace", font_size=22, color=GRAY_B)
            .next_to(row, DOWN, buff=0.18)
            for ch, row in zip(message, bit_rows)
        ])
        self.play(Write(bit_labels))
        self.next_slide()

        all_bits = bytes_to_bits(msg_bytes)
        full_row = bit_row(all_bits, size=0.44)
        full_row.move_to(ORIGIN)

        # ℓ → \ell
        ell_lbl = MathTex(
            r"\ell = %d \text{ bits (message length)}" % len(all_bits),
            font_size=28, color=WHITE
        )
        ell_lbl.next_to(full_row, DOWN, buff=0.5)

        self.play(
            FadeOut(char_group), FadeOut(bit_labels),
            *[ReplacementTransform(r.copy(), full_row) for r in bit_rows],
            FadeOut(bit_rows),
        )
        self.play(Write(ell_lbl))
        self.next_slide()

        # FIXED: Tex → MathTex for math content
        note = MathTex(
            r"\text{SHA-256 can hash messages up to } 2^{64}-1 \text{ bits long}",
            font_size=22, color=GRAY_B
        )
        note.next_to(ell_lbl, DOWN, buff=0.4)
        self.play(FadeIn(note, shift=UP * 0.1))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 2 — Preprocessing / Padding
# ============================================================================

class Preprocessing(Slide):
    def construct(self):
        message   = "abc"
        msg_bytes = message.encode("ascii")
        ell       = len(msg_bytes) * 8  # 24
        k         = 423
        TW        = 12.2   # target visual width

        title = section_title("Step 2 · Padding to 512 bits", color=TEAL)
        self.play(Write(title))

        msg_bits = bytes_to_bits(msg_bytes)

        row1 = bit_row(msg_bits, size=0.44)
        row1.move_to(UP * 1.6)
        cap1 = MathTex(
            r"\text{``abc''} \rightarrow \ell = 24 \text{ bits}",
            font_size=26, color=GRAY_B
        )
        cap1.next_to(row1, UP, buff=0.28)
        self.play(FadeIn(row1), Write(cap1))
        self.next_slide()

        def make_full_row(msg, one=True, zeros=True, length=True, scale_to=TW):
            cells = VGroup(*[bit_sq(b, size=0.44) for b in msg])
            if one:
                ob = bit_sq('1', size=0.44)
                ob.set_stroke(color=AMBER, width=3)
                cells.add(ob)
            cells.arrange(RIGHT, buff=0.04)
            g = VGroup(cells)
            if zeros:
                zb = Rectangle(width=0.44 * 8, height=0.44,
                               fill_color=BIT_OFF, fill_opacity=1,
                               stroke_color="#505068", stroke_width=0.5)
                zb.next_to(cells, RIGHT, buff=0.05)
                g.add(zb)
            if length:
                lb_bits = int_to_bits(ell, 64)
                lb = VGroup(*[bit_sq(b, size=0.44, on="#9D4EDD") for b in lb_bits])
                lb.arrange(RIGHT, buff=0.04)
                lb.next_to(g[-1], RIGHT, buff=0.05)
                g.add(lb)
            if scale_to:
                g.set(width=scale_to)
            g.move_to(UP * 1.6)
            return g

        row2 = make_full_row(msg_bits, one=True,  zeros=False, length=False)
        cap2 = MathTex(r"\text{Append single ``1'' bit}", font_size=26, color=AMBER)
        cap2.next_to(row2, UP, buff=0.28)
        self.play(FadeOut(row1), FadeOut(cap1))
        self.play(FadeIn(row2), Write(cap2))
        self.next_slide()

        row3 = make_full_row(msg_bits, one=True,  zeros=True,  length=False)
        # FIXED: Tex → MathTex, \pmod now in math mode
        cap3 = MathTex(
            rf"\text{{Append }} k = {k} \text{{ zero bits  (so }} \ell+1+k \equiv 448 \pmod{{512}} \text{{)}}",
            font_size=24, color=GRAY_B
        )
        cap3.next_to(row3, UP, buff=0.28)
        self.play(FadeOut(row2), FadeOut(cap2))
        self.play(FadeIn(row3), Write(cap3))
        self.next_slide()

        row4 = make_full_row(msg_bits, one=True,  zeros=True,  length=True)
        # FIXED: Tex → MathTex, \ell is math
        cap4 = MathTex(
            rf"\text{{Append }} \ell = {ell} \text{{ as 64-bit big-endian integer}}",
            font_size=24, color="#9D4EDD"
        )
        cap4.next_to(row4, UP, buff=0.28)
        self.play(FadeOut(row3), FadeOut(cap3))
        self.play(FadeIn(row4), Write(cap4))
        self.next_slide()

        total = MathTex(
            r"\text{Total} = 24 + 1 + 423 + 64 = 512 \text{ bits} = \text{one 512-bit block } M^{(1)}",
            font_size=26, color=AMBER
        )
        total.next_to(row4, DOWN, buff=0.55)
        self.play(Write(total))
        self.next_slide()

        self.play(FadeOut(VGroup(row4, cap4, total)))

        padded = pad_message(msg_bytes)
        block  = to_blocks(padded)[0]
        words  = [int.from_bytes(block[i * 4:i * 4 + 4], 'big') for i in range(16)]

        word_boxes = VGroup()
        for i, w in enumerate(words):
            box_group = VGroup(
                hex_box(w, width=1.64, height=0.64, font_size=21),
                MathTex(f"M_{{{i}}}", font_size=17, color=GRAY_B),
            ).arrange(DOWN, buff=0.07)
            word_boxes.add(box_group)
        word_boxes.arrange_in_grid(rows=4, cols=4, buff=0.22)
        word_boxes.move_to(DOWN * 0.25)

        wcap = MathTex(
            r"M^{(1)} \text{ parsed into sixteen 32-bit words } M_0 \ldots M_{15}",
            font_size=25, color=SKY
        )
        wcap.next_to(word_boxes, UP, buff=0.42)

        self.play(Write(wcap))
        self.play(LaggedStart(*[FadeIn(w, scale=0.85) for w in word_boxes], lag_ratio=0.07))

        # FIXED: Tex → MathTex, \ell is math
        note = MathTex(
            r"M_0 = \mathtt{0x61626380} \text{ encodes 'abc' + the 1-bit} \qquad "
            r"M_{15} = \mathtt{0x00000018} = 24 = \ell",
            font_size=22, color=GRAY_B
        )
        note.to_edge(DOWN, buff=0.45)
        self.play(FadeIn(note, shift=UP * 0.1))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 3 — IV and Constants: where they come from
# ============================================================================

class IVAndConstants(Slide):
    def construct(self):
        title = section_title("Step 3 · Initial Values & Round Constants", color=AMBER)
        self.play(Write(title))

        iv_sub = MathTex(
            r"\text{Initial Hash Values } H_0^{(0)} \ldots H_7^{(0)}",
            font_size=30, color=AMBER
        )
        iv_sub.next_to(title, DOWN, buff=0.45)
        self.play(Write(iv_sub))
        self.next_slide()

        primes8 = [2, 3, 5, 7, 11, 13, 17, 19]
        import math

        derive_rows = VGroup()
        for i, (p, h) in enumerate(zip(primes8, H0_INIT)):
            sqrt_p = math.sqrt(p)
            p_txt = MathTex(r"\sqrt{%d}" % p, font_size=22, color=TEAL)
            eq = Text("→  frac part → first 32 bits →", font_size=20, color=GRAY_B)
            h_txt = MathTex(
                r"H_{%d}^{(0)} = \mathtt{%08x}" % (i, h),
                font_size=22, color=AMBER
            )
            row = VGroup(p_txt, eq, h_txt).arrange(RIGHT, buff=0.35)
            derive_rows.add(row)

        derive_rows.arrange(DOWN, buff=0.20, aligned_edge=LEFT)
        derive_rows.next_to(iv_sub, DOWN, buff=0.50)
        derive_rows.shift(LEFT * 0.8)

        self.play(LaggedStart(
            *[FadeIn(r, shift=RIGHT * 0.25) for r in derive_rows],
            lag_ratio=0.10))
        self.next_slide()

        insight_iv = Text(
            "Nothing-up-my-sleeve numbers: the choice of H₀…H₇ is publicly verifiable\n"
            "and provably not hiding any backdoor.",
            font_size=22, color=GRAY_B, line_spacing=1.4, slant=ITALIC)
        insight_iv.to_edge(DOWN, buff=0.55)
        self.play(FadeIn(insight_iv, shift=UP * 0.15))
        self.next_slide()

        self.play(FadeOut(VGroup(iv_sub, derive_rows, insight_iv)))

        k_sub = MathTex(
            r"\text{Round Constants } K_0 \ldots K_{63}",
            font_size=30, color=LAVENDER
        )
        k_sub.next_to(title, DOWN, buff=0.45)
        self.play(Write(k_sub))

        primes64 = []
        candidate = 2
        while len(primes64) < 64:
            if all(candidate % p for p in range(2, int(candidate**0.5) + 1)):
                primes64.append(candidate)
            candidate += 1

        # FIXED: Tex → MathTex, contains \sqrt[3]
        kderive = MathTex(
            r"K_t = \text{first 32 bits of the fractional part of } \sqrt[3]{p_t}\\ "
            r"\text{where } p_0=2, p_1=3, p_2=5, \ldots p_{63}=311",
            font_size=24, color=GRAY_B
        )
        kderive.next_to(k_sub, DOWN, buff=0.55)
        self.play(FadeIn(kderive, shift=UP * 0.15))
        self.next_slide()

        k_boxes = VGroup()
        for i in range(8):
            bx = hex_box(K_CONSTANTS[i], width=1.58, height=0.72,
                         font_size=20, color=LAVENDER)
            lbl = MathTex(f"K_{{{i}}}", font_size=16, color=GRAY_B)
            grp = VGroup(lbl, bx).arrange(DOWN, buff=0.06)
            k_boxes.add(grp)
        k_boxes.arrange(RIGHT, buff=0.22)
        k_boxes.next_to(kderive, DOWN, buff=0.50)

        kdots = MathTex(
            r"\ldots K_8 \text{ through } K_{63} \text{ derived the same way } \ldots",
            font_size=21, color=GRAY_B
        )
        kdots.next_to(k_boxes, DOWN, buff=0.30)

        self.play(LaggedStart(*[FadeIn(bx, scale=0.85) for bx in k_boxes], lag_ratio=0.07))
        self.play(Write(kdots))

        insight_k = Text(
            "64 unique constants ensure every round is cryptographically distinct.",
            font_size=22, color=AMBER, slant=ITALIC)
        insight_k.to_edge(DOWN, buff=0.55)
        self.play(FadeIn(insight_k, shift=UP * 0.15))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 4 — Merkle–Damgård construction
# ============================================================================

class MerkleDamgard(Slide):
    def construct(self):
        title = section_title("Step 4 · Merkle–Damgård Construction", color=LAVENDER)
        self.play(Write(title))

        sub = Text(
            "How SHA-256 handles arbitrarily long messages one 512-bit block at a time",
            font_size=24, color=GRAY_B)
        sub.next_to(title, DOWN, buff=0.35)
        self.play(FadeIn(sub, shift=UP * 0.15))
        self.next_slide()

        BOX_W, BOX_H = 2.0, 1.1
        Y_CHAIN      = 0.6

        def compress_box(label, x, col=SKY):
            bx = RoundedRectangle(corner_radius=0.12, width=BOX_W, height=BOX_H,
                                   fill_color=DKGRAY, fill_opacity=1,
                                   stroke_color=col, stroke_width=2.5)
            lt = MathTex(label, font_size=22, color=col)
            lt.move_to(bx.get_center())
            grp = VGroup(bx, lt)
            grp.move_to(RIGHT * x + UP * Y_CHAIN)
            return grp

        iv_box = RoundedRectangle(corner_radius=0.10, width=1.3, height=BOX_H,
                                   fill_color=DKGRAY, fill_opacity=1,
                                   stroke_color=AMBER, stroke_width=2)
        iv_lbl = MathTex(r"\text{IV}\\ H^{(0)}", font_size=20, color=AMBER)
        iv_lbl.move_to(iv_box.get_center())
        iv_grp = VGroup(iv_box, iv_lbl).move_to(LEFT * 5.5 + UP * Y_CHAIN)

        c1 = compress_box(r"\text{compress}\\ M^{(1)}", -2.8, SKY)
        c2 = compress_box(r"\text{compress}\\ M^{(2)}",  0.4, SKY)
        c3 = compress_box(r"\text{compress}\\ M^{(N)}",  3.6, SKY)

        dots = Text("⋯", font_size=36, color=GRAY_B).move_to(RIGHT * 2.0 + UP * Y_CHAIN)

        digest_box = RoundedRectangle(corner_radius=0.10, width=1.5, height=BOX_H,
                                       fill_color=DKGRAY, fill_opacity=1,
                                       stroke_color=TEAL, stroke_width=2)
        digest_lbl = MathTex(r"\text{digest}\\ H^{(N)}", font_size=20, color=TEAL)
        digest_lbl.move_to(digest_box.get_center())
        digest_grp = VGroup(digest_box, digest_lbl).next_to(c3, RIGHT, buff=0.05)

        def block_label(txt, ref):
            t = Text(txt, font_size=19, color=GRAY_B)
            t.next_to(ref, DOWN, buff=0.65)
            return t

        blk1_lbl = block_label("512-bit\nblock 1", c1)
        blk2_lbl = block_label("512-bit\nblock 2", c2)
        blkN_lbl = block_label("512-bit\nblock N", c3)

        def h_arrow(src, dst, col=WHITE):
            return Arrow(src.get_right(), dst.get_left(),
                         buff=0.12, color=col, stroke_width=2.5, tip_length=0.18)

        arr0 = h_arrow(iv_grp, c1, AMBER)
        arr1 = h_arrow(c1, c2, SKY)
        arr2 = h_arrow(c2, dots, SKY)
        arr3 = h_arrow(dots, c3, SKY)

        def feed_arrow(dst_box, src_txt, col=SKY):
            return Arrow(src_txt.get_top() + UP * 0.05, dst_box.get_bottom(),
                         buff=0.0, color=col, stroke_width=2.5, tip_length=0.15)

        fa1 = feed_arrow(c1, blk1_lbl, SKY)
        fa2 = feed_arrow(c2, blk2_lbl, SKY)
        fa3 = feed_arrow(c3, blkN_lbl, SKY)

        self.play(FadeIn(iv_grp, scale=0.9))
        self.play(GrowArrow(arr0), FadeIn(c1, scale=0.9))
        self.play(FadeIn(fa1, shift=UP * 0.2), FadeIn(blk1_lbl))
        self.next_slide()

        self.play(GrowArrow(arr1), FadeIn(c2, scale=0.9))
        self.play(FadeIn(fa2, shift=UP * 0.2), FadeIn(blk2_lbl))
        self.next_slide()

        self.play(GrowArrow(arr2), FadeIn(dots))
        self.play(GrowArrow(arr3), FadeIn(c3, scale=0.9))
        self.play(FadeIn(fa3, shift=UP * 0.2), FadeIn(blkN_lbl))
        self.play(FadeIn(digest_grp, shift=LEFT * 0.2))
        self.next_slide()

        le_note = Text(
            "The chained structure means every output bit depends on\n"
            "every preceding block — a change anywhere cascades forward.",
            font_size=20, color=GRAY_B, line_spacing=1.3)
        le_note.to_edge(DOWN, buff=1.6)
        self.play(FadeIn(le_note, shift=UP * 0.15))
        self.next_slide()

        # FIXED: Tex → MathTex, contains H^{(i)} and math
        dm_note = MathTex(
            r"\text{Each 'compress' step uses a Davies--Meyer feedforward:}\\ "
            r"H^{(i)} = \text{compress}(H^{(i-1)}, M^{(i)}) + H^{(i-1)} \quad "
            r"\text{(more on this soon)}",
            font_size=19, color=LAVENDER
        )
        dm_note.next_to(le_note, DOWN, buff=0.25)
        self.play(FadeIn(dm_note, shift=UP * 0.1))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 5 — Message Schedule  W₀ … W₆₃
# ============================================================================

class MessageSchedule(Slide):
    def construct(self):
        _, _, _, blocks, _, all_W, _ = sha256_full(b"abc")
        W = all_W[0]

        # Title with LaTeX
        title_tex = MathTex(r"\text{Step 5 · Message Schedule } W_0 \ldots W_{63}", font_size=36, color=SKY)
        title_tex.to_edge(UP, buff=0.40)
        self.play(Write(title_tex))
        sub = MathTex(
            r"\text{The 16 message words are expanded to 64 using } \sigma_0 \text{ and } \sigma_1",
            font_size=25, color=GRAY_B
        )
        sub.next_to(title_tex, DOWN, buff=0.35)
        self.play(FadeIn(sub, shift=UP * 0.15))
        self.next_slide()

        ROW_H, COL_W, COLS = 0.70, 1.55, 8
        start_y = 1.3

        def make_w_box(i, color=SKY, bg=DKGRAY, sw=1.8):
            bx  = RoundedRectangle(corner_radius=0.07, width=COL_W, height=ROW_H,
                                    fill_color=bg, fill_opacity=1,
                                    stroke_color=color, stroke_width=sw)
            idx = MathTex(f"W_{{{i}}}", font_size=16, color=GRAY_B)
            val = Text(f"{W[i]:08x}", font="monospace", font_size=18)
            idx.move_to(bx.get_top()    + DOWN * 0.18)
            val.move_to(bx.get_center() + DOWN * 0.10)
            grp = VGroup(bx, idx, val)
            col = i % COLS
            row = i // COLS
            grp.move_to(
                RIGHT * (col - (COLS - 1) / 2) * (COL_W + 0.17) +
                UP    * (start_y - row * (ROW_H + 0.20)))
            return grp

        boxes_16 = VGroup(*[make_w_box(i) for i in range(16)])
        w16_cap = MathTex(
            r"W_0 \ldots W_{15} \text{ (taken directly from the 512-bit message block)}",
            font_size=21, color=SKY
        )
        w16_cap.next_to(boxes_16, UP, buff=0.25)

        self.play(Write(w16_cap))
        self.play(LaggedStart(*[FadeIn(bx, scale=0.85) for bx in boxes_16],
                               lag_ratio=0.05), run_time=1.8)
        self.next_slide()

        formula = MathTex(
            r"W_t = \sigma_1(W_{t-2}) + W_{t-7} + \sigma_0(W_{t-15}) + W_{t-16}"
            r"\qquad (16 \le t \le 63)",
            font_size=28, color=AMBER)
        formula.next_to(boxes_16, DOWN, buff=0.45)
        self.play(Write(formula))

        sigma_defs = MathTex(
            r"\sigma_0(x) = \mathrm{ROTR}^7(x) \oplus \mathrm{ROTR}^{18}(x) \oplus \mathrm{SHR}^3(x)",
            font_size=22, color=TEAL)
        sigma1_def = MathTex(
            r"\sigma_1(x) = \mathrm{ROTR}^{17}(x) \oplus \mathrm{ROTR}^{19}(x) \oplus \mathrm{SHR}^{10}(x)",
            font_size=22, color=CORAL)
        VGroup(sigma_defs, sigma1_def).arrange(DOWN, buff=0.18)
        sigma_grp = VGroup(sigma_defs, sigma1_def)
        sigma_grp.next_to(formula, DOWN, buff=0.28)
        self.play(FadeIn(sigma_grp, shift=DOWN * 0.15))
        self.next_slide()

        sources = {0: CORAL, 1: TEAL, 14: SKY, 15: LAVENDER}
        self.play(*[boxes_16[idx][0].animate.set_stroke(col, width=4)
                    for idx, col in sources.items()])

        w16_box = make_w_box(16, color=AMBER, bg="#2A1A4A", sw=3)
        w16_box.next_to(formula, DOWN, buff=0.35).shift(LEFT * 4.5)
        self.play(FadeOut(sigma_grp), FadeIn(w16_box, scale=1.2))

        def curved_arr(src, dst, col):
            return CurvedArrow(src.get_bottom(), dst.get_top(),
                               angle=-0.5, color=col,
                               stroke_width=2.5, tip_length=0.14)

        arrows = VGroup(
            curved_arr(boxes_16[15], w16_box, LAVENDER),
            curved_arr(boxes_16[14], w16_box, SKY),
            curved_arr(boxes_16[1],  w16_box, TEAL),
            curved_arr(boxes_16[0],  w16_box, CORAL),
        )
        self.play(LaggedStart(*[Create(a) for a in arrows], lag_ratio=0.2))
        self.next_slide()

        extra = VGroup()
        for i in range(17, 20):
            bx = make_w_box(i, color=AMBER, bg="#2A1A4A", sw=2)
            bx.next_to(w16_box, RIGHT, buff=0.18 + (i - 17) * (COL_W + 0.18))
            extra.add(bx)
        self.play(LaggedStart(*[FadeIn(bx, scale=0.85) for bx in extra], lag_ratio=0.1))

        ff_lbl = MathTex(r"\ldots \text{ up to } W_{63}", font_size=22, color=GRAY_B)
        ff_lbl.next_to(extra, RIGHT, buff=0.3)
        w63_box = make_w_box(63, color=CORAL, bg="#2A1A4A", sw=3)
        w63_box.next_to(ff_lbl, RIGHT, buff=0.3)
        self.play(Write(ff_lbl), FadeIn(w63_box, scale=1.2))
        self.next_slide()

        key = MathTex(
            r"\text{Every round gets a distinct } W_t \text{ that incorporates ALL 512 bits of the block}",
            font_size=22, color=AMBER
        )
        key.to_edge(DOWN, buff=0.45)
        self.play(Write(key))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 6 — Bitwise operations: Ch, Maj, Σ, σ
# ============================================================================

class BitwiseOperations(Slide):
    def construct(self):
        _, _, _, _, all_rounds, _, _ = sha256_full(b"abc")
        r0      = all_rounds[0][0]
        e, f, g = r0['before'][4], r0['before'][5], r0['before'][6]
        a, b, c = r0['before'][0], r0['before'][1], r0['before'][2]
        ch_val  = Ch(e, f, g)
        maj_val = Maj(a, b, c)

        DEMO = 10
        eb   = format(e,       '032b')[-DEMO:]
        fb   = format(f,       '032b')[-DEMO:]
        gb   = format(g,       '032b')[-DEMO:]
        chb  = format(ch_val,  '032b')[-DEMO:]
        ab   = format(a,       '032b')[-DEMO:]
        bb_s = format(b,       '032b')[-DEMO:]
        cb   = format(c,       '032b')[-DEMO:]
        majb = format(maj_val, '032b')[-DEMO:]

        SQ, VSEP = 0.40, 0.72

        title = section_title("Step 6 · Bitwise Logic Inside SHA-256", color=CORAL)
        sub = MathTex(
            r"\text{Ch, Maj, } \Sigma_0, \Sigma_1, \sigma_0, \sigma_1 \text{ --- the six mixing functions}",
            font_size=25, color=GRAY_B
        )
        sub.next_to(title, DOWN, buff=0.32)
        self.play(Write(title))
        self.play(FadeIn(sub, shift=UP * 0.15))
        self.next_slide()

        # ── Ch ──
        ch_head = VGroup(
            MathTex(
                r"\text{Ch}(e,f,g) = (e \land f) \oplus (\lnot e \land g)",
                font_size=28, color=AMBER
            ),
            Text('"Choose": where e=1 output f, where e=0 output g',
                 font_size=22, color=GRAY_B),
        ).arrange(DOWN, buff=0.14)
        ch_head.next_to(sub, DOWN, buff=0.50)
        self.play(Write(ch_head[0]), FadeIn(ch_head[1], shift=UP * 0.1))
        self.next_slide()
        anchor = ch_head.get_bottom() + DOWN * 0.55

        OFF = "#202035"
        def bit_row_group(label, bits, color, y_off):
            lbl = Text(label, font_size=26, color=color, font="monospace")
            row = brow(bits, sq=SQ, on=color, off=OFF)
            grp = VGroup(lbl, row).arrange(RIGHT, buff=0.4)
            grp.move_to(anchor + DOWN * y_off)
            return grp, row

        ge, re = bit_row_group("e", eb, CORAL,  0)
        gf, rf = bit_row_group("f", fb, TEAL,   VSEP)
        gg, rg = bit_row_group("g", gb, SKY,    VSEP * 2)
        for grp in (ge, gf, gg):
            self.play(FadeIn(grp, shift=RIGHT * 0.25), run_time=0.5)
        self.next_slide()

        rule = Text("Column rule: Ch[i] = f[i] if e[i]=1, else g[i]",
                    font_size=21, color=AMBER)
        rule.next_to(gg, DOWN, buff=0.42)
        self.play(Write(rule))

        ch_cells = VGroup()
        for i in range(DEMO):
            cell = Square(side_length=SQ,
                          fill_color=((AMBER if fb[i] == "1" else OFF) if eb[i] == "1" else (AMBER if gb[i]=="1" else OFF)),
                          fill_opacity=1, stroke_color="#505068", stroke_width=0.5)
            ch_cells.add(cell)
        ch_cells.arrange(RIGHT, buff=0.05)
        ch_lbl = Text("Ch", font_size=24, color=AMBER, font="monospace")
        ch_row_grp = VGroup(ch_lbl, ch_cells).arrange(RIGHT, buff=0.4)
        ch_row_grp.next_to(gg, DOWN, buff=VSEP - 0.05).align_to(ge, LEFT)
        ch_cells.align_to(re, LEFT)

        col_anims = [
            AnimationGroup(
                Indicate(rf[i] if eb[i] == "1" else rg[i], color=AMBER, scale_factor=1.35, run_time=0.3),
                FadeIn(ch_cells[i], run_time=0.3),
            ) for i in range(DEMO)
        ]
        self.play(FadeOut(rule))
        self.play(Write(ch_lbl), LaggedStart(*col_anims, lag_ratio=0.15))
        ch_hex = MathTex(
            r"= \mathtt{0x%08x} \text{ (low %d bits: } %s\text{)}" % (ch_val, DEMO, chb),
            font_size=21, color=AMBER
        )
        ch_hex.next_to(ch_row_grp, DOWN, buff=0.28)
        self.play(Write(ch_hex))
        self.next_slide()

        self.play(FadeOut(VGroup(ch_head, ge, gf, gg, ch_row_grp, ch_hex)))

        # ── Maj ──
        maj_head = VGroup(
            MathTex(
                r"\text{Maj}(a,b,c) = (a \land b) \oplus (a \land c) \oplus (b \land c)",
                font_size=28, color=LAVENDER
            ),
            Text('"Majority vote": output 1 when ≥2 of the 3 inputs are 1',
                 font_size=22, color=GRAY_B),
        ).arrange(DOWN, buff=0.14)
        maj_head.next_to(sub, DOWN, buff=0.50)
        self.play(Write(maj_head[0]), FadeIn(maj_head[1], shift=UP * 0.1))

        ga, ra   = bit_row_group("a", ab,   CORAL,  0)
        gb2, rb2 = bit_row_group("b", bb_s, TEAL,   VSEP)
        gc, rc   = bit_row_group("c", cb,   SKY,    VSEP * 2)
        for grp in (ga, gb2, gc):
            self.play(FadeIn(grp, shift=RIGHT * 0.25), run_time=0.5)
        self.next_slide()

        maj_cells = VGroup()
        for i in range(DEMO):
            ones = int(ab[i]) + int(bb_s[i]) + int(cb[i])
            cell = Square(side_length=SQ,
                          fill_color=(LAVENDER if ones >= 2 else "#202035"),
                          fill_opacity=1, stroke_color="#505068", stroke_width=0.5)
            maj_cells.add(cell)
        maj_cells.arrange(RIGHT, buff=0.05)
        maj_lbl = Text("Maj", font_size=24, color=LAVENDER, font="monospace")
        maj_row_grp = VGroup(maj_lbl, maj_cells).arrange(RIGHT, buff=0.4)
        maj_row_grp.next_to(gc, DOWN, buff=VSEP - 0.05).align_to(ga, LEFT)
        maj_cells.align_to(ra, LEFT)

        vote_anims = [
            AnimationGroup(
                *[Indicate(x, color=LAVENDER, scale_factor=1.3, run_time=0.3)
                  for x, bit in zip([ra[i], rb2[i], rc[i]], [ab[i], bb_s[i], cb[i]])
                  if bit == "1"],
                FadeIn(maj_cells[i], run_time=0.3),
            ) for i in range(DEMO)
        ]
        self.play(Write(maj_lbl), LaggedStart(*vote_anims, lag_ratio=0.15))
        maj_hex = MathTex(
            r"= \mathtt{0x%08x} \text{ (low %d bits: } %s\text{)}" % (maj_val, DEMO, majb),
            font_size=21, color=LAVENDER
        )
        maj_hex.next_to(maj_row_grp, DOWN, buff=0.28)
        self.play(Write(maj_hex))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])

        # ── Sigma / sigma rotation demo ──
        rot_title = MathTex(
            r"\text{Rotation Functions } \Sigma_0, \Sigma_1, \sigma_0, \sigma_1",
            font_size=38, color=WHITE
        )
        rot_title.to_edge(UP, buff=0.40)
        rot_note = Text(
            "Each diffuses bits by XOR-ing three rotated (or shifted) copies of one word.",
            font_size=25, color=GRAY_B
        )
        rot_note.next_to(rot_title, DOWN, buff=0.32)
        self.play(Write(rot_title))
        self.play(FadeIn(rot_note, shift=UP * 0.15))
        self.next_slide()

        table_data = [
            (r"\Sigma_0(a)", r"\mathrm{ROTR}^2(a) \oplus \mathrm{ROTR}^{13}(a) \oplus \mathrm{ROTR}^{22}(a)", AMBER),
            (r"\Sigma_1(e)", r"\mathrm{ROTR}^6(e) \oplus \mathrm{ROTR}^{11}(e) \oplus \mathrm{ROTR}^{25}(e)", CORAL),
            (r"\sigma_0(x)", r"\mathrm{ROTR}^7(x) \oplus \mathrm{ROTR}^{18}(x) \oplus \mathrm{SHR}^3(x)", TEAL),
            (r"\sigma_1(x)", r"\mathrm{ROTR}^{17}(x) \oplus \mathrm{ROTR}^{19}(x) \oplus \mathrm{SHR}^{10}(x)", SKY),
        ]
        table_rows = VGroup()
        for name, defn, col in table_data:
            nt = MathTex(name, font_size=24, color=col)
            eq = Text("=", font_size=22, color=WHITE)
            dt = MathTex(defn, font_size=22, color=GRAY_B)
            row = VGroup(nt, eq, dt).arrange(RIGHT, buff=0.35)
            table_rows.add(row)
        table_rows.arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        table_rows.next_to(rot_note, DOWN, buff=0.55)
        self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.2) for r in table_rows],
                               lag_ratio=0.15))
        self.next_slide()

        demo_val  = e & 0xFFFF
        demo_bits = format(demo_val, '016b')

        demo_lbl  = Text("e  (low 16 bits)", font_size=23, color=CORAL)
        demo_row  = brow(demo_bits, sq=0.37, on=CORAL, off="#202035")
        demo_grp  = VGroup(demo_lbl, demo_row).arrange(RIGHT, buff=0.4)
        demo_grp.to_edge(DOWN, buff=2.2)
        self.play(FadeIn(demo_grp, shift=DOWN * 0.2))
        self.next_slide()

        n6       = 6
        rot6     = ((demo_val >> n6) | (demo_val << (16 - n6))) & 0xFFFF
        rot6_bits = format(rot6, '016b')

        rot_lbl  = MathTex(r"\mathrm{ROTR}^6(e)", font_size=23, color=AMBER)
        rot_row  = brow(rot6_bits, sq=0.37, on=AMBER, off="#202035")
        rot_grp  = VGroup(rot_lbl, rot_row).arrange(RIGHT, buff=0.4)
        rot_grp.next_to(demo_grp, UP, buff=0.45).align_to(demo_grp, LEFT)

        moving = VGroup(*[cell.copy().set_color(AMBER).set_fill(AMBER, 1)
                          for cell in demo_row])
        self.add(moving)
        self.play(*[
            moving[i].animate.move_to(rot_row[(i + n6) % 16].get_center())
            for i in range(16)
        ], run_time=1.5)
        self.remove(moving)
        self.play(FadeIn(rot_grp))

        explain = Text(
            f"rotate right {n6} positions → rightmost {n6} bits wrap to the front",
            font_size=20, color=GRAY_B)
        explain.next_to(demo_grp, DOWN, buff=0.28)
        self.play(Write(explain))
        self.next_slide()

        insight = Text(
            "Three rotations chosen so every output bit is a function of all input bits —\n"
            "no bit position is left 'invisible' to any other.",
            font_size=21, color=AMBER, line_spacing=1.4, slant=ITALIC)
        insight.next_to(explain, DOWN, buff=0.5)
        self.play(FadeOut(explain), FadeIn(insight, shift=UP * 0.15))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 7 — The 64-round compression function
# ============================================================================

class CompressionRound(Slide):
    def construct(self):
        digest, H_final, padded, blocks, all_rounds, all_W, all_H_in = sha256_full(b"abc")
        rounds  = all_rounds[0]
        H_in_0  = all_H_in[0]

        TITLE_Y  =  3.70
        RECIPE_Y =  1.85
        REG_Y    = -0.20
        PANEL_Y  = -1.70
        INFO_Y   = -3.25

        title = Text("Step 7 · SHA-256 Compression — 64 rounds", font_size=36, color=AMBER)
        title.move_to(UP * TITLE_Y)
        self.play(Write(title))

        f_t1 = MathTex(
            r"T_1 = h + \Sigma_1(e) + \mathrm{Ch}(e,f,g) + K_t + W_t",
            font_size=27, color=T1_COLOR)
        f_t2 = MathTex(
            r"T_2 = \Sigma_0(a) + \mathrm{Maj}(a,b,c)",
            font_size=27, color=T2_COLOR)
        f_new = MathTex(
            r"a' = T_1 + T_2 \qquad e' = d + T_1",
            font_size=27, color=NEW_COLOR)
        f_shift = MathTex(
            r"(b',c',d',f',g',h') \;\leftarrow\; (a,b,c,e,f,g)",
            font_size=24, color=GRAY_B)
        recipe = VGroup(f_t1, f_t2, f_new, f_shift).arrange(DOWN, buff=0.14)
        recipe.move_to(UP * RECIPE_Y)
        self.play(Write(f_t1), Write(f_t2))
        self.play(Write(f_new), Write(f_shift))
        self.next_slide()

        regs = make_reg_row(H_in_0, REG_Y)
        self.play(FadeIn(regs))

        def info_text(t, r):
            return Text(
                f"Round t={t:2d}   W{t} = {r['W']:08x}   K{t} = {r['K']:08x}",
                font="monospace", font_size=22,
            ).move_to(UP * INFO_Y)

        info = info_text(0, rounds[0])
        self.play(Write(info))
        self.next_slide()

        def do_round(t, narrate=True):
            r       = rounds[t]
            new_info = info_text(t, r)
            anims   = [Transform(info, new_info)]
            if narrate:
                anims += [
                    regs[7].box.animate.set_stroke(T1_COLOR, width=4),
                    regs[4].box.animate.set_stroke(T1_COLOR, width=4),
                    regs[5].box.animate.set_stroke(T1_COLOR, width=4),
                    regs[6].box.animate.set_stroke(T1_COLOR, width=4),
                    regs[0].box.animate.set_stroke(T2_COLOR, width=4),
                    regs[1].box.animate.set_stroke(T2_COLOR, width=4),
                    regs[2].box.animate.set_stroke(T2_COLOR, width=4),
                ]
            self.play(*anims, run_time=(0.9 if narrate else 0.08))

            def temp_box(label, val, col, pos):
                bx  = RoundedRectangle(corner_radius=0.08, width=3.5, height=0.82,
                                        stroke_color=col, stroke_width=3,
                                        fill_color=DKGRAY, fill_opacity=1)
                row = VGroup(
                    MathTex(label + " =", font_size=22, color=col),
                    Text(f"{val:08x}", font="monospace", font_size=24),
                ).arrange(RIGHT, buff=0.2)
                row.move_to(bx.get_center())
                grp = VGroup(bx, row)
                grp.move_to(pos)
                return grp

            positions = [LEFT * 3.0 + UP * PANEL_Y, RIGHT * 3.0 + UP * PANEL_Y]
            t1_grp = temp_box("T_1", r['T1'], T1_COLOR, positions[0])
            t2_grp = temp_box("T_2", r['T2'], T2_COLOR, positions[1])

            if narrate:
                self.play(FadeIn(t1_grp, shift=UP * 0.25), FadeIn(t2_grp, shift=UP * 0.25))
            else:
                self.add(t1_grp, t2_grp)

            new_vals = [None] * 8
            shift_a  = []
            for src, dst in [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7)]:
                mv = regs[src].val_t.copy()
                shift_a.append(mv.animate.move_to(regs[dst].val_t.get_center()))
                new_vals[dst] = mv

            shift_a.append(regs[3].val_t.animate.move_to(t1_grp.get_center()).set_opacity(0))
            shift_a.append(regs[7].val_t.animate.move_to(t1_grp.get_center()).set_opacity(0))
            shift_a.append(regs[0].val_t.animate.move_to(t2_grp.get_center()).set_opacity(0))

            new_a = Text(f"{r['after'][0]:08x}", font="monospace",
                          font_size=REG_FONT, color=NEW_COLOR)
            new_a.move_to(t2_grp.get_center())
            new_e = Text(f"{r['after'][4]:08x}", font="monospace",
                          font_size=REG_FONT, color=NEW_COLOR)
            new_e.move_to(t1_grp.get_center())
            shift_a.append(new_a.animate.move_to(regs[0].val_t.get_center()))
            shift_a.append(new_e.animate.move_to(regs[4].val_t.get_center()))
            new_vals[0] = new_a;  new_vals[4] = new_e

            self.add(new_a, new_e)
            self.play(*shift_a, run_time=1.3 if narrate else 0.10)

            for i in range(8):
                self.remove(regs[i].val_t)
            self.remove(new_a, new_e)
            for i in range(8):
                nv = new_vals[i]
                nv.set_opacity(1)
                if i not in (0, 4):
                    nv.set_color(WHITE)
                regs[i].val_t         = nv
                regs[i].submobjects[2] = nv
                self.add(nv)

            self.play(FadeOut(t1_grp), FadeOut(t2_grp),
                      run_time=(0.8 if narrate else 0.06))
            if narrate:
                self.play(*[regs[i].box.animate.set_stroke(WHITE, width=2)
                             for i in range(8)])
                self.play(*[new_vals[i].animate.set_color(WHITE) for i in (0, 4)])

        do_round(0, narrate=True);  self.next_slide()
        do_round(1, narrate=True);  self.next_slide()

        ff = Text("… 62 more rounds (t = 2 … 63), same recipe each time …",
                   font_size=26, color=GRAY_B).move_to(UP * PANEL_Y)
        self.play(FadeIn(ff))
        
        final_state = rounds[63]['after']
        jump_anims = []
        
        for i in range(8):
            new_val = Text(f"{final_state[i]:08x}", font="monospace",
                           font_size=REG_FONT, color=WHITE)
            new_val.move_to(regs[i].val_t.get_center())
            jump_anims.append(Transform(regs[i].val_t, new_val))
            regs[i].val_t = new_val
            regs[i].submobjects[2] = new_val

        self.play(*jump_anims, run_time=1.0)

        self.play(FadeOut(ff), FadeOut(info))
        self.next_slide()

        self.play(FadeOut(recipe))
        h0_row     = make_reg_row(H_in_0, RECIPE_Y)
        h0_caption = MathTex(
            r"H^{(0)} \text{ (initial hash value for this block)}",
            font_size=24, color=GRAY_B
        )
        h0_caption.next_to(h0_row, UP, buff=0.14)
        plus = Text("+", font_size=36).move_to((h0_row.get_center() + regs.get_center()) / 2)
        self.play(FadeIn(h0_row), Write(h0_caption), Write(plus))
        self.next_slide()

        h1_row     = make_reg_row(H_final, INFO_Y + 1.4)
        h1_caption = MathTex(r"H^{(1)} = \text{final digest!}", font_size=24, color=AMBER)
        h1_caption.next_to(h1_row, DOWN, buff=0.24)
        eq = Text("=", font_size=36).move_to((regs.get_center() + h1_row.get_center()) / 2)
        self.play(FadeIn(h1_row), Write(h1_caption), Write(eq))
        self.next_slide()

        digest_line = VGroup(
            Text('SHA-256("abc")  =', font_size=21, color=GRAY_B),
            Text(digest.hex(), font="monospace", font_size=19, color=AMBER),
        ).arrange(DOWN, buff=0.10)
        digest_line.to_edge(DOWN, buff=0.30)
        self.play(Write(digest_line))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 8 — Davies–Meyer construction (why it's one-way)
# ============================================================================

class DaviesMeyer(Slide):
    def construct(self):
        title = section_title(
            "Step 8 · Davies–Meyer: \nThe Feedforward That Makes It One-Way",
            color=TEAL, size=34
        )
        self.play(Write(title))
        self.next_slide()

        eq = MathTex(
            r"H^{(i)} = \underbrace{E_{M^{(i)}}(H^{(i-1)})}_{\text{block cipher run as compressor}}"
            r"\;+\; H^{(i-1)}",
            font_size=34, color=WHITE
        )
        eq.next_to(title, DOWN, buff=0.50)
        self.play(Write(eq))
        self.next_slide()

        y_main = 0.4

        hin_box = RoundedRectangle(
            corner_radius=0.10, width=1.5, height=0.75,
            fill_color=DKGRAY, fill_opacity=1,
            stroke_color=AMBER, stroke_width=2
        )
        hin_lbl = MathTex(r"H^{(i-1)}", font_size=20, color=AMBER)
        hin_lbl.move_to(hin_box.get_center())
        hin_grp = VGroup(hin_box, hin_lbl).move_to(LEFT * 5.0 + UP * y_main)

        em_box = RoundedRectangle(
            corner_radius=0.10, width=2.5, height=0.75,
            fill_color="#2A2A3C", fill_opacity=1,
            stroke_color=SKY, stroke_width=2.5
        )
        em_lbl = VGroup(
            MathTex(r"E_{M^{(i)}}", font_size=22, color=SKY),
            Text("(64 SHA-256 rounds)", font_size=15, color=GRAY_B),
        ).arrange(DOWN, buff=0.05)
        em_lbl.move_to(em_box.get_center())
        em_grp = VGroup(em_box, em_lbl).move_to(LEFT * 1.5 + UP * y_main)

        xor_circ = Circle(
            radius=0.34, color=LAVENDER, stroke_width=2.5,
            fill_color=DKGRAY, fill_opacity=1
        )
        xor_txt = Text("⊕", font_size=26, color=LAVENDER)
        xor_txt.move_to(xor_circ.get_center())
        xor_grp = VGroup(xor_circ, xor_txt).move_to(RIGHT * 2.5 + UP * y_main)

        hout_box = RoundedRectangle(
            corner_radius=0.10, width=1.5, height=0.75,
            fill_color=DKGRAY, fill_opacity=1,
            stroke_color=TEAL, stroke_width=2
        )
        hout_lbl = MathTex(r"H^{(i)}", font_size=20, color=TEAL)
        hout_lbl.move_to(hout_box.get_center())
        hout_grp = VGroup(hout_box, hout_lbl).move_to(RIGHT * 4.8 + UP * y_main)

        msg_lbl = MathTex(r"M^{(i)} \text{ (used as block-cipher key)}", font_size=18, color=GRAY_B)
        msg_lbl.next_to(em_grp, DOWN, buff=0.30)
        msg_arr = Arrow(
            msg_lbl.get_top(), em_grp.get_bottom(),
            buff=0.05, color=GRAY_B, stroke_width=2, tip_length=0.13
        )

        ff_path = CurvedArrow(
            hin_grp.get_bottom() + DOWN * 0.05,
            xor_grp.get_bottom() + DOWN * 0.05,
            angle=0.85,
            color=AMBER, stroke_width=2.5, tip_length=0.15
        )

        arr_main1 = Arrow(
            hin_grp.get_right(), em_grp.get_left(),
            buff=0.1, color=WHITE, stroke_width=2.5, tip_length=0.15
        )
        arr_main2 = Arrow(
            em_grp.get_right(), xor_grp.get_left(),
            buff=0.1, color=WHITE, stroke_width=2.5, tip_length=0.15
        )
        arr_main3 = Arrow(
            xor_grp.get_right(), hout_grp.get_left(),
            buff=0.1, color=TEAL, stroke_width=2.5, tip_length=0.15
        )

        self.play(FadeIn(hin_grp, scale=0.9), FadeIn(em_grp, scale=0.9))
        self.play(FadeIn(msg_lbl), GrowArrow(msg_arr))
        self.play(GrowArrow(arr_main1), GrowArrow(arr_main2))
        self.play(FadeIn(xor_grp, scale=0.9))
        self.play(Create(ff_path))
        self.play(GrowArrow(arr_main3), FadeIn(hout_grp, scale=0.9))
        self.next_slide()

        diagram_mobjects = VGroup(
            hin_grp, em_grp, xor_grp, hout_grp,
            msg_lbl, msg_arr,
            ff_path, arr_main1, arr_main2, arr_main3
        )
        self.play(FadeOut(diagram_mobjects))

        ow_bullets = VGroup(
            Text(
                "▸  Without the feedforward (+H_in), reversing the block cipher\n"
                "    would directly invert the hash — game over.",
                font_size=21, color=GRAY_B, line_spacing=1.3
            ),
            Text(
                "▸  With the feedforward, even a perfectly invertible cipher\n"
                "    leaves the attacker needing to solve x + E(x) = target\n"
                "    — computationally infeasible.",
                font_size=21, color=GRAY_B, line_spacing=1.3
            ),
            Text(
                "▸  This is the Davies–Meyer mode, proven secure in the\n"
                "    ideal-cipher model (Shannon, 1949 / Davies-Meyer, 1985).",
                font_size=21, color=LAVENDER, line_spacing=1.3
            ),
        )
        ow_bullets.arrange(DOWN, buff=0.42, aligned_edge=LEFT)
        ow_bullets.next_to(eq, DOWN, buff=0.55)
        ow_bullets.shift(LEFT * 0.3)

        for bul in ow_bullets:
            self.play(FadeIn(bul, shift=RIGHT * 0.2))
            self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 9 — Avalanche effect
# ============================================================================

class AvalancheEffect(Slide):
    def construct(self):
        msg1 = b"abc"
        msg2 = bytes([msg1[0] ^ 0x01]) + msg1[1:]

        d1, _, _, _, _, _, _ = sha256_full(msg1)
        d2, _, _, _, _, _, _ = sha256_full(msg2)

        bits1     = format(int.from_bytes(d1, 'big'), '0256b')
        bits2     = format(int.from_bytes(d2, 'big'), '0256b')
        diff_mask = [b1 != b2 for b1, b2 in zip(bits1, bits2)]
        n_diff    = sum(diff_mask)

        title = section_title("Step 9 · The Avalanche Effect", color=CORAL)
        title.to_edge(UP, buff=0.35)
        sub = subtitle_text("Change 1 input bit → roughly half of all output bits flip")
        sub.next_to(title, DOWN, buff=0.18)
        self.play(Write(title))
        self.play(FadeIn(sub, shift=UP * 0.15))
        self.next_slide()

        m1l = MathTex(r"\text{msg}_1 = \text{``abc''}", font_size=28, color=TEAL)
        m2l = MathTex(r"\text{msg}_2 = \text{``\`bc''}", font_size=28, color=CORAL)
        dn  = Text("(bit 0 of 'a' flipped: 0x61 → 0x60)",
                    font_size=20, color=GRAY_B)
        VGroup(m1l, m2l, dn).arrange(DOWN, buff=0.20).next_to(sub, DOWN, buff=0.45)
        self.play(Write(m1l), Write(m2l))
        self.play(FadeIn(dn, shift=UP * 0.1))
        self.next_slide()

        bits_a1  = byte_to_bits(ord('a'))
        bits_a2  = byte_to_bits(ord('`'))
        SQ2      = 0.44

        def msg_bit_row(bits, col):
            cells = VGroup(*[
                Square(side_length=SQ2,
                       fill_color=(col if b == "1" else "#202035"),
                       fill_opacity=1, stroke_color="#505060", stroke_width=0.8)
                for b in bits
            ])
            cells.arrange(RIGHT, buff=0.06)
            return cells

        r1 = msg_bit_row(bits_a1, TEAL)
        r2 = msg_bit_row(bits_a2, CORAL)
        VGroup(r1, r2).arrange(DOWN, buff=0.22)

        flip_idx = next(i for i, (x, y) in enumerate(zip(bits_a1, bits_a2)) if x != y)
        l1 = Text("'a' = 0x61", font_size=20, color=TEAL)
        l2 = Text("'`' = 0x60", font_size=20, color=CORAL)
        VGroup(l1, l2).arrange(DOWN, buff=0.22)
        bit_display = VGroup(VGroup(l1, l2), VGroup(r1, r2)).arrange(RIGHT, buff=0.35)
        bit_display.next_to(dn, DOWN, buff=0.45)
        self.play(FadeIn(bit_display, shift=DOWN * 0.2))
        self.next_slide()

        flip1, flip2 = r1[flip_idx], r2[flip_idx]
        for _ in range(3):
            self.play(
                flip1.animate.set_stroke(AMBER, width=5),
                flip2.animate.set_stroke(AMBER, width=5), run_time=0.25)
            self.play(
                flip1.animate.set_stroke("#505060", width=0.8),
                flip2.animate.set_stroke("#505060", width=0.8), run_time=0.20)
        self.next_slide()

        self.play(FadeOut(bit_display), FadeOut(dn),
                  FadeOut(m1l), FadeOut(m2l))

        SQD  = 0.245
        BUFD = 0.028
        COLS = 32
        BLOCK_CY = -0.35
        ROW_H    = 8 * (SQD + BUFD)
        LBL_H    = 0.45
        GAP      = 0.55
        y_top =  BLOCK_CY + GAP / 2 + ROW_H / 2
        y_bot =  BLOCK_CY - GAP / 2 - ROW_H / 2

        def digest_grid(bits, diff, label, label_col, y_centre):
            cells = VGroup()
            for i, b in enumerate(bits):
                clr = (CORAL if b == "1" else "#5A1020") if diff[i] \
                      else (TEAL  if b == "1" else "#202035")
                stk = CORAL if diff[i] else "#404055"
                cells.add(Square(side_length=SQD, fill_color=clr, fill_opacity=1,
                                 stroke_color=stk, stroke_width=0.5))
            cells.arrange_in_grid(cols=COLS, buff=BUFD)
            cells.move_to(UP * y_centre)
            lbl = Text(label, font_size=18, color=label_col)
            lbl.next_to(cells, UP, buff=0.14)
            return VGroup(lbl, cells), cells

        grp1, g1 = digest_grid(bits1, [False] * 256,
                                f'SHA-256("abc")  =  {d1.hex()[:16]}…', TEAL,  y_top)
        grp2, g2 = digest_grid(bits2, diff_mask,
                                f'SHA-256("`bc")  =  {d2.hex()[:16]}…', CORAL, y_bot)

        self.play(LaggedStart(*[FadeIn(c, scale=0.7) for c in g1], lag_ratio=0.006),
                  run_time=1.5)
        self.play(Write(grp1[0]))
        self.next_slide()

        self.play(
            LaggedStart(*[
                AnimationGroup(FadeIn(c, scale=1.4, run_time=0.4)) if diff_mask[i]
                else FadeIn(c, run_time=0.25)
                for i, c in enumerate(g2)
            ], lag_ratio=0.005),
            run_time=2.0)
        self.play(Write(grp2[0]))
        self.next_slide()

        diff_info = Text(
            f"{n_diff} / 256 bits differ  ({n_diff / 256 * 100:.1f}%)",
            font_size=24, color=AMBER)
        diff_info.to_edge(DOWN, buff=0.40)

        legend = VGroup(
            Square(side_length=0.22, fill_color=TEAL,  fill_opacity=1, stroke_width=0),
            Text("same",      font_size=17, color=TEAL),
            Square(side_length=0.22, fill_color=CORAL, fill_opacity=1, stroke_width=0),
            Text("different", font_size=17, color=CORAL),
        ).arrange(RIGHT, buff=0.16)
        legend.next_to(diff_info, UP, buff=0.18)

        self.play(FadeIn(legend), Write(diff_info))
        diff_cells = [g2[i] for i in range(256) if diff_mask[i]]
        self.play(*[c.animate.set_stroke(AMBER, width=3) for c in diff_cells], run_time=0.7)
        self.play(*[c.animate.set_stroke(CORAL, width=0.5) for c in diff_cells], run_time=0.5)
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# SCENE 10 — Final Summary / Security properties
# ============================================================================

class FinalSummary(Slide):
    def construct(self):
        title = section_title("SHA-256: Why It's Secure", color=WHITE)
        self.play(Write(title))

        # All descriptions now use MathTex with \text{} for text
        props = [
            ("One-way (preimage resistance)",
             r"\text{Given digest } d \text{, finding } m \text{ with } \text{SHA-256}(m) = d \text{ requires } \sim 2^{256} \text{ operations}",
             AMBER),
            ("Second preimage resistance",
             r"\text{Given } m_1 \text{, finding } m_2 \neq m_1 \text{ with the same digest requires } \sim 2^{256} \text{ work}",
             TEAL),
            ("Collision resistance",
             r"\text{Finding any } m_1 \neq m_2 \text{ with the same digest requires } \sim 2^{128} \text{ work (birthday)}",
             CORAL),
            ("Avalanche",
             r"\text{Every output bit depends on every input bit -- 1-bit change flips } \sim 50\% \text{ of digest}",
             LAVENDER),
            ("Deterministic",
             r"\text{Same input always gives the same digest -- no hidden randomness}",
             SKY),
            ("Fixed output size",
             r"\text{Any message length from } 0 \text{ to } 2^{64}-1 \text{ bits } \rightarrow \text{ always exactly 256 bits out}",
             AMBER),
        ]

        prop_group = VGroup()
        for name, desc, clr in props:
            hdr = Text(f"▸  {name}", font_size=25, color=clr)
            bod = MathTex(desc, font_size=19, color=GRAY_B)
            row = VGroup(hdr, bod).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
            prop_group.add(row)

        prop_group.arrange(DOWN, buff=0.42, aligned_edge=LEFT)
        prop_group.next_to(title, DOWN, buff=0.60)
        prop_group.shift(LEFT * 0.4)

        for row in prop_group:
            self.play(FadeIn(row[0], shift=RIGHT * 0.3), run_time=0.45)
            self.play(FadeIn(row[1], shift=RIGHT * 0.2), run_time=0.35)
            self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])

        closing = VGroup(
            Text("SHA-256", font_size=80, color=AMBER, weight=BOLD),
            Text("Used in TLS, Bitcoin, Git, Docker, code-signing, and more.",
                 font_size=26, color=GRAY_B),
            Text("Based on FIPS PUB 180-4.",
                 font_size=22, color=GRAY_B, slant=ITALIC),
        ).arrange(DOWN, buff=0.55)
        closing.move_to(ORIGIN)
        self.play(Write(closing[0], run_time=1.2))
        self.play(FadeIn(closing[1], shift=UP * 0.15))
        self.play(FadeIn(closing[2], shift=UP * 0.10))
        self.next_slide()

        self.play(*[FadeOut(m) for m in self.mobjects])


# ============================================================================
# COMBINED PRESENTATION
# ============================================================================

class SHA256Presentation(
    TitleSlide,
    MessageToBits,
    Preprocessing,
    IVAndConstants,
    MerkleDamgard,
    MessageSchedule,
    BitwiseOperations,
    CompressionRound,
    DaviesMeyer,
    AvalancheEffect,
    FinalSummary,
):
    """
    Combines all scenes into one continuous manim-slides presentation.
    """
    def construct(self):
        TitleSlide.construct(self)
        MessageToBits.construct(self)
        Preprocessing.construct(self)
        IVAndConstants.construct(self)
        MerkleDamgard.construct(self)
        MessageSchedule.construct(self)
        BitwiseOperations.construct(self)
        CompressionRound.construct(self)
        DaviesMeyer.construct(self)
        AvalancheEffect.construct(self)
        FinalSummary.construct(self)
