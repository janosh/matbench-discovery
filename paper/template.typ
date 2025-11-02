#let arkheion(
  title: "",
  abstract: [],
  keywords: (),
  authors: (),
  date: none,
  body,
) = {
  // basic document properties
  set document(
    author: authors.map(author => if type(author.name) == str {
      return author.name
    } else if author.name.has("text") {
      return author.name.text
    } else {
      return ""
    }),
    title: title,
  )
  set page(
    margin: (left: 22mm, right: 22mm, top: 22mm, bottom: 25mm),
    numbering: "1",
    number-align: center,
  )
  set text(font: "New Computer Modern", lang: "en")
  show math.equation: set block(spacing: 0.65em)
  set math.equation(numbering: "(1)")
  set heading(numbering: "1.1")

  // Set run-in subheadings, starting at level 4.
  show heading: it => {
    // H1 and H2
    if it.level == 1 {
      pad(
        bottom: 10pt,
        it,
      )
    } else if it.level == 2 {
      pad(
        bottom: 8pt,
        it,
      )
    } else if it.level > 3 {
      text(11pt, weight: "bold", it.body + " ")
    } else {
      it
    }
  }

  // Title
  line(length: 100%, stroke: 2pt)
  pad(
    bottom: 4pt,
    top: 4pt,
    align(center)[
      #block(text(weight: 500, 1.75em, title))
      #v(1em, weak: true)
    ],
  )
  line(length: 100%, stroke: 2pt)

  // Author info
  v(2em)
  pad(x: 1em)[
    #grid(
      columns: (1fr,) * calc.min(5, authors.len()),
      row-gutter: 2em,
      column-gutter: 1em,
      ..authors.map(author => align(center)[
        #if author.keys().contains("orcid") {
          link("http://orcid.org/" + author.orcid)[
            #grid(
              columns: (auto, 1pt),
              text(size: 0.8em)[*#author.name*], [#pad(left: 6pt, top: -4pt, image("figs/orcid.svg", width: 8pt))],
            )
          ]
        } else { [*#author.name*] }
        #v(-6pt)
        #if author.keys().contains("email") { text(size: 0.7em, weight: 100)[#author.email] }
        #if author.keys().contains("affiliation") { text(size: 0.7em, weight: 100)[#author.affiliation] }
      ]),
    )]

  v(2em)
  align(center)[#date]

  // Abstract.
  pad(
    x: 2em,
    top: 1em,
    bottom: 0.4em,
    align(center)[
      #heading(
        outlined: false,
        numbering: none,
        text(0.85em, smallcaps[Abstract]),
      )
      #set par(justify: true)
      #set text(hyphenate: false)

      #abstract
    ],
  )

  // Keywords
  if keywords.len() > 0 {
    [*_Keywords_* #h(0.3cm)] + keywords.map(str).join(" · ")
  }

  // Main body
  // equations: reference as "eq. (1)"
  set math.equation(numbering: "(1)", supplement: none)
  show ref: it => {
    // wrap equation numbers in parentheses when referencing
    if it.element != none and it.element.func() == math.equation {
      link(it.target)[eq.~(#it)]
    } else {
      it
    }
  }

  // dark blue links and references
  show ref: set text(fill: blue.darken(20%))
  show link: set text(fill: blue.darken(20%))

  // change sub/superscript font size
  set sub(size: 0.8em)
  set super(size: 0.8em)

  // style tables
  set table(
    inset: (x: 5pt, y: 4pt), // cell padding
    // blue shade for header row, light gray for first column
    fill: (col, row) => if row == 0 { blue.lighten(90%) } else if col == 0 { luma(245) } else { none },
    // thin horizontal lines between rows (except header), none between columns
    stroke: (_, y) => if y > 0 { (top: 0.2pt) },
  )
  // bold table headers
  show table.cell.where(y: 0): set text(weight: "bold")

  // paragraphs
  set par(justify: true, first-line-indent: 1em)

  // figures
  show figure: set text(size: 0.95em)

  body
}

#let arkheion-appendices(body) = {
  counter(heading).update(0)
  counter("appendices").update(1)

  set heading(
    numbering: (..nums) => {
      let vals = nums.pos()
      let value = "ABCDEFGHIJ".at(vals.at(0) - 1)
      if vals.len() == 1 {
        return "APPENDIX " + value
      } else {
        return value + "." + nums.pos().slice(1).map(str).join(".")
      }
    },
  )
  [#pagebreak() #body]
}

#let subfigure-kind = "subfigure"

#let subfigure-counter = counter(subfigure-kind)

#let subfigure(
  body,
  pos: bottom + center,
  dx: 0%,
  dy: 6%,
  caption: "",
  numbering: "a)",
  separator: none,
  label: none,
  supplement: none,
  placement: top,
) = {
  subfigure-counter.step()

  let fig = figure(
    body,
    caption: none,
    kind: subfigure-kind,
    supplement: none,
    numbering: numbering,
    outlined: false,
    placement: placement,
  )

  if caption != "" and separator == none {
    separator = ":"
  }

  context {
    let sub-fig-num = subfigure-counter.display(numbering)
    let caption-content = [#supplement #sub-fig-num#separator #caption]
    return [ #fig#label #place(pos, dx: dx, dy: dy, caption-content) ]
  }
}

// reset subfigure counter when out of the parent figure
#show figure: itm => {
  if itm.kind != subfigure-kind {
    subfigure-counter.update(0)
  }
  itm
}

// Custom rule for formatting references to subfigures as "Figure 1a)"
#show ref: itm => {
  let elem = itm.element
  // Check if referenced element is a subfigure
  if elem != none and elem.func() == figure and elem.kind == subfigure-kind {
    // Find all outlined figures before the subfigure's location
    let outlined-figs-before = query(figure.where(outlined: true).before(elem.location()))
    // Filter out figures that are tables (kind: table)
    let actual-figs-before = outlined-figs-before.filter(f => f.kind != table)
    // The parent figure is the last one in the filtered list
    let parent-fig = actual-figs-before.last()
    // Calculate the parent figure number based on its position in the filtered list
    let parent-num = actual-figs-before.len()
    // Get the subfigure counter state array at the element's location
    let subfig-state = subfigure-counter.at(elem.location())
    // Format the state using the numbering function
    let subfig-num = numbering(elem.numbering, ..subfig-state)
    // Note: This assumes standard '1, 2, 3, ...' numbering for parent figures.
    // Custom parent numbering formats (e.g., "A.1") won't be replicated here.
    return [#parent-fig.supplement #parent-num#subfig-num]
  }
  // Default handling for all other references
  itm
}
