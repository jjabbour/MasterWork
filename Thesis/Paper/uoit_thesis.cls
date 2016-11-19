\ProvidesClass{uoit_thesis}[2007/04/16 UOIT Thesis Class]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}

\ProcessOptions

\LoadClass[letterpaper]{report}
%\LoadClass[letterpaper,openright,12pt,titlepage,final]{book}
\long\def\comment#1{}
\comment{

  Example of use:
    \documentclass[letterpaper,12pt,MSc]{uoit_thesis}
    \begin{document}
    \title{How to Write Thesis\\
            With Two Line Titles}
    \author{John Henry Candidate}
    \principaladviser{John Parker}
    \firstreader{John Green}
    \secondreader{John BigBooty\\(Another Department)}

    \beforeabstract
    \prefacesection{Abstract}
        This thesis answers all the questions you were afraid to ask ...
    \afterabstract
    \prefacesection{Acknowledgements}
        I would like to thank...
    \afterpreface

    \chapter{Introduction}
         ...
    \chapter{Conclusions}
         ...
    \appendix
    \chapter{A Long Proof}
         ...
%    \bibliographystyle{alpha}
%    \bibliography{mybib}


Documentation:
    It sets the margins and interline spacing and
    disallows page breaks at hyphens.

    The \beforeabstract command creates the title page, a copyright
    page (default), and a signature page (optional), the tables of
    contents, tables and figures. The signature page is not mentioned
    in the rules for submission, and by implication are not allowed.

    The user then includes the preface using the
    \prefacesection{Abstract} command.

    Next should come an \afterabstract command to insert the declaration.

    Then the user should put preface section(s), using the
    \prefacesection{section title} command. The first of these should
    be the statement of education and research, others may be
    acknowledgements etc.

    The \afterpreface command then  sets things up to start
    the main body (continuing the sequence of page numbering, but
    putting page headings etc.).

    Having one of \beforepreface or \afterpreface without the other
    may cause problems, so don't do it.

    The following commands can control what goes in the front matter
    material:

        \documentclass[letterpaper,12pt]{uoit_thesis}
        \title{thesis title}
        \author{author's name}
        \dept{author's department}
                - Computer Science if omitted
        \principaladviser{the principal advisor's name}
           (or \principaladvisor, if you prefer advisor spelled with o)
        \firstreader{the first reader's name}
        \secondreader{the second reader's name}
        \submitdate{month year in which submitted to GPO}
                - date LaTeX'd if omitted
        \copyrightyear{year degree conferred (next year if submitted in Dec.)}
                - year LaTeX'd (or next year, in December) if omitted
        \copyrighttrue or \copyrightfalse
                - produce or don't produce a copyright page (true by default)
        \figurespagetrue or \figurespagefalse
                - produce or don't produce a List of Figures page
                  (true by default)
        \tablespagetrue or \tablespagefalse
                - produce or don't produce a List of Tables page
                  (true by default)
        \sigpagetrue or \sigpagefalse
                - produce or don't produce a signature page
                  (again false by default, as Manchester theses don't
                   require this. However, it's left in so that you can
                   your own `examiners page' if you like!)

        If you set all these to false, then perhaps you don't want any
        preface at all, so perhaps you hould leave out both \beforepreface
        and \afterpreface -).

This style uses interline spacing that is 1.5 times normal (about double
spaced), except in the figure and table environments where normal spacing
is used.

That can be changed by doing:
    \renewcommand{\baselinestretch}{1.6}
(or whatever you want instead of 1.5)
Unfortunately, that command probably won't take effect unless it comes
before the \begin{document} or after something like \small \normalsize.

The example above shows the 12pt substyle being used.  This seems to give
acceptable looking results, but it may be omitted to get smaller print.

} %end of comment


% We need >15mm  margins except on the binding edge, where it is >40mm
%\oddsidemargin 14.6truemm \evensidemargin 0truemm
%\marginparwidth 40pt \marginparsep 10pt
%\topmargin 0truemm \headsep 10truemm
%\textheight 224truemm \textwidth 137truemm

%%%
%  the setup for the printer maps2
%%%
\setlength{\voffset}{-2mm}
\setlength{\oddsidemargin}{12.0mm}
\setlength{\evensidemargin}{12.0mm}
%-.4 mm
\setlength{\textwidth}{154mm}
\setlength{\topmargin}{-11mm}
% 0mm !!!
\setlength{\textheight}{228mm}

%%%
%  the setup for the printer eps3
%%%
%\setlength{\voffset}{-2mm}
%\setlength{\oddsidemargin}{17.6mm}
%\setlength{\evensidemargin}{-.4 mm}
%\setlength{\textwidth}{140mm}
%\setlength{\topmargin}{1.5mm}
%\setlength{\textheight}{228mm}

%%%
%  the setup for printer macol
%%%
%\setlength{\voffset}{-2mm}
%\setlength{\oddsidemargin}{17.6mm}
%\setlength{\evensidemargin}{-.4 mm}
%\setlength{\textwidth}{140mm}
%\setlength{\topmargin}{3mm}
%\setlength{\textheight}{228mm}


% Disallow page breaks at hyphens (this will give some underfull vbox's,
% so an alternative is to use \brokenpenalty=100 and manually search
% for and fix such page breaks)
\brokenpenalty=100

% Use 1.5 times the normal baseline-to-baseline skip (approx double space)
%%%%%%%%%%%%%%%%%%%%%%%  this is the line spacing%%%%%%%%%%%%%%%%%%%%
\renewcommand{\baselinestretch}{1.5}



% Set the headings style and numbering of equations
\usepackage{fancyhdr,amsmath}
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]%
               {\markboth{\chaptername\ \thechapter.\ #1}{}}
\renewcommand{\sectionmark}[1]%
             {\markright{ \thesection.\ #1}}
\fancyhf{}
\fancyhead[LO]{\nouppercase{\slshape \rightmark}}
\fancyhead[RE]{\nouppercase{\slshape \leftmark}}
\cfoot{\thepage}

\numberwithin{equation}{section}
\renewcommand{\theequation}%
{\arabic{chapter}.\arabic{section}.\arabic{equation}}

% Redefine the macro used for floats (including figures and tables)
% so that single spacing is used.
% (Note \def\figure{\@float{figure}set single spacing} doesn't work
%  because figure has an optional argument)
\def\@xfloat#1[#2]{\ifhmode \@bsphack\@floatpenalty -\@Mii\else
   \@floatpenalty-\@Miii\fi\def\@captype{#1}\ifinner
      \@parmoderr\@floatpenalty\z@
    \else\@next\@currbox\@freelist{\@tempcnta\csname ftype@#1\endcsname
       \multiply\@tempcnta\@xxxii\advance\@tempcnta\sixt@@n
       \@tfor \@tempa :=#2\do
                        {\if\@tempa h\advance\@tempcnta \@ne\fi
                         \if\@tempa t\advance\@tempcnta \tw@\fi
                         \if\@tempa b\advance\@tempcnta 4\relax\fi
                         \if\@tempa p\advance\@tempcnta 8\relax\fi
         }\global\count\@currbox\@tempcnta}\@fltovf\fi
    \global\setbox\@currbox\vbox\bgroup
    \def\baselinestretch{1}\@normalsize
    \boxmaxdepth\z@
    \hsize\columnwidth \@parboxrestore}

% Redefine the macro used for footnotes to use single spacing
\long\def\@footnotetext#1{\insert\footins{\def\baselinestretch{1}\footnotesize
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
   \edef\@currentlabel{\csname p@footnote\endcsname\@thefnmark}\@makefntext
    {\rule{\z@}{\footnotesep}\ignorespaces
      #1\strut}}}

% \author, \title are defined in report; here are the rest of the
% front matter defining macros
% Note: If you really want, you can define your own macro for, say
% \supervisor , just by copying the definitions of \principaladvisor
%
\def\dept#1{\gdef\@dept{#1}}
\def\faculty#1{\gdef\@faculty{#1}}
\def\degreetitle#1{\gdef\@degreetitle{#1}}
\def\submitiondate#1{\gdef\@submitiondate{#1}}
\def\wordcount#1{\gdef\@wordcount{#1}}
% \author, \title and \date defined in report
\def\@title{}\def\@author{}\def\@dept{Faculty of Science, University of Ontario Institute of Technology}
\def\@faculty{Faculty of Science}\def\@degreetitle{Master of Science}
\def\@date{\ifcase\the\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space \number\the\year}
\def\@submitiondate{\ifcase\the\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space \number\the\year}

\newif\iffigurespage \newif\iftablespage
\figurespagetrue \tablespagetrue

\def\makepreliminaries{\pagenumbering{roman}
                       \pagestyle{plain}}

%% Creates the titlepage
%\def\titlepage{%
%        \thispagestyle{empty}%
%        %\null\vskip1in%
%        \begin{center}
%                \Huge\uppercase\expandafter{\@title}
%        \end{center}
%        \vfill
%        \begin{figure*}[h!]
%           \centering
%           \includegraphics[width=1.4in]{uoit_logo.eps}
%        \end{figure*}
%        \vfill
%        \begin{center}
%                \rm By\\
%                \@author\\
%                \@dept
%        \end{center}
%                \vfill
%        \begin{center}
%                \@submitiondate
%        \end{center}
%        \vfill
%        \begin{center}
%        \sc A thesis submitted to the \\
%        University of Ontario Institute of Technology \\
%        in accordance with the requirements of the degree\\
%        of \expandafter{\@degreetitle}  in the \expandafter{\@faculty}
%        \end{center}
%
%	\vspace{1.5cm}
%	\begin{flushright}
%	%\@wordcount
%	\end{flushright}
%
%        \newpage}

% Creates the titlepage
\def\titlepage{%
        \thispagestyle{empty}%
        %\null\vskip1in%
        \begin{center}
                \Huge\uppercase\expandafter{\@title}
        \end{center}
        \vfill
%        \begin{center}
%           \epsfig{file=uoit_logo.eps, width=1.8in}
%        \end{center}
%        \vfill
        \begin{center}
                \rm By\\
                \@author\\
        \end{center}
                \vfill
        \begin{center}

        \end{center}
        \vfill
        \begin{center}
		A Thesis Submitted in Partial Fulfillment\\ of the Requirements for the Degree of \\
        \expandafter{\@degreetitle} \\  in \\ The \expandafter{\@faculty} \\
        \end{center}
        \vfill
        \begin{center}
         Modelling and Computational Science
        \end{center}
        \vfill
   	    \begin{center}
        August 2015
        \end{center}
        \begin{center}
        \textcopyright \@author, 2015
        \end{center}


	\vspace{1.5cm}
	\begin{flushright}
	%\@wordcount
	\end{flushright}

        \newpage}








% Creates the page with the author's declaration
\def\declarationpage{%
\prefacesection{Author's Declaration}
        \begin{center}
        \begin{minipage}{4in}
           \parindent=0pt I declare that this work
           was carried out in accordance with the regulations of
           the University of Ontario \linebreak Institute of Technology. The
	   work is original except where indicated by special reference
	   in the text and no part of this document has been
	   submitted for any other degree.  Any views expressed in
	   the dissertation are those of the author and in no way
	   represent those of the \linebreak University of Ontario Institute
	   of Technology.  This document has not been presented to
	   any other University for examination either in Canada
	   or overseas.



        \vspace{1in}

        \hbox to 4in{\hfil\shortstack{\vrule width 3in height
            0.4pt\\ \expandafter{\@author}}}
        \vspace{1in}
        \parindent=1in Date: \expandafter{July $4$, $2015$}
%\@date
                \par
        \end{minipage}%
        \end{center}%
        \vfill\newpage}

%Creates the table of contents with list of figures and tables as options
% The \addvspace business below is to avoid the extra space that gets put
% in the tables of figures and tables between chapters (it makes no sense
% now that the numbering is continuous)

\def\contents{\pagestyle{fancy}
             \tableofcontents
              \newpage
              %  \iftablespage
              %   \addvspace{10pt}
              %    \listoftables
              %   \newpage
              % \fi
               \iffigurespage
                \renewcommand{\baselinestretch}{0.6}
                \addvspace{10pt}
                \listoffigures
                \listoftables
                \newpage
        \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%NEU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates Krupp acknowledgement
\def\krupppage{%
\prefacesection{Financial Support}
        \begin{center}
        \begin{minipage}{4.3in}
           \parindent=0pt

        \vspace{2.0cm}

        {\Large Diese Arbeit...}

        \vspace{3.0cm}

        This dissertation...

        \end{minipage}%
        \end{center}%
        \vfill\newpage}
%%%%%%%%%%%%%%%%%%%%%%%%%%NEU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\afterpreliminaries{\pagenumbering{arabic}
               \pagestyle{fancy}}



\def\prefacesection#1{%
        \chapter*{#1}
        \addcontentsline{toc}{chapter}{#1}}


%%% DPC 15/9/93 Make the bibliography and the index generate an entry
%%% in the table of contents.
%\renewcommand\appendix{\par
%  \setcounter{chapter}{0}%
%  \setcounter{section}{0}%
%  \renewcommand\@chapapp{\appendixname}%
%  \renewcommand\thechapter{\@Alph\c@chapter}
%  \renewcommand\chaptername{Appendix}
%
%\fancyhf{}
%\fancyhead[LO]{\nouppercase{\slshape \leftmark}}
%\fancyhead[RE]{\nouppercase{\slshape \leftmark}}
%\cfoot{\thepage}
%
%
% \renewcommand{\theequation}{\Alph{chapter}.\arabic{equation}}}
%
%\newdimen\bibindent
%\setlength\bibindent{1.5em}
%\renewenvironment{thebibliography}[1]
%     {\chapter*{\bibname
%        \@mkboth{\bibname}{\bibname}}%
%             \addcontentsline{toc}{chapter}{Bibliography}%
%             \renewcommand{\baselinestretch}{1}\@normalsize%
%         \list{\@biblabel{\@arabic\c@enumiv}}%
%           {\settowidth\labelwidth{\@biblabel{#1}}%
%            \leftmargin\labelwidth
%            \advance\leftmargin\labelsep
%            \@openbib@code
%            \usecounter{enumiv}%
%            \let\p@enumiv\@empty
%            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
%      \sloppy
%      \clubpenalty4000
%      \@clubpenalty \clubpenalty
%      \widowpenalty4000%
%      \sfcode`\.\@m}
%     {\def\@noitemerr
%       {\@latex@warning{Empty `thebibliography' environment}}%
%      \endlist}

\renewenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\@makeschapterhead{\indexname}]%
                \@mkboth{\indexname}{\indexname}%
                \addcontentsline{toc}{chapter}{Index}%
                \renewcommand{\baselinestretch}{1}\@normalsize%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
