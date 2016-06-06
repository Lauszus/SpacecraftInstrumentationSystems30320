# Spacecraft Instrumentation Systems 30320
#### June, 2016
_________
[![Build Status](https://travis-ci.com/Lauszus/SpacecraftInstrumentationSystems30320.svg?token=ppc6rHRAs23cjxNFyjc1&branch=master)](https://travis-ci.com/Lauszus/SpacecraftInstrumentationSystems30320)

If you use Git, then please do __NOT__ commit to the master branch. Instead make a separate branch and send a [pull request](https://help.github.com/articles/using-pull-requests), then it will be merged into master if there is no conflicts. I would recommend the following GUI if you are new to Git: <https://www.sourcetreeapp.com>.

The ShareLaTeX project is located here: <https://www.sharelatex.com/project/57308fe73f0401db1321aae2>.

An overview of the different subjects can be found at the following [Google spreadsheet](https://docs.google.com/spreadsheets/d/1nXihh6wFuoOudE7xgI5fpLbBgiGJ6B4esyN-JFeMGO0). Please note that people who has done presentations on the respective subjects, gets first priority on these for the report. If you see that someone is working on the same as you, please contact him/her and see if you can merge your materials together, so we do not end up with duplicates in the report.

If you are new to LaTex, then I would recommend the following editor: <http://www.xm1math.net/texmaker>. You will also need a Latex distribution: [Windows](http://miktex.org), [Mac](https://tug.org/mactex) and [Linux](http://www.tug.org/texlive).

Please also join the Facebook group: <https://www.facebook.com/groups/1711333799143909>.

Our YouTube channels is located at the following link: <https://www.youtube.com/channel/UC-DMik5iodpwCzGbywo2CMg>. If you need to add a video, then please contact me and I will upload it.

### The book

The final version of the book can be downloaded at the following link: [Spacecraft Instrumentation Systems - Europa Life Finder Mission.pdf](https://github.com/Lauszus/SpacecraftInstrumentationSystems30320/releases/download/1.0.0/Spacecraft.Instrumentation.Systems.-.Europa.Life.Finder.Mission.pdf).

### Notes

Please use proper citations for references instead of footnotes. More information can be found at the following [link](https://www.sharelatex.com/learn/Bibliography_management_with_natbib).

Also see the examples that is already added to the book: [bibliography/biblio.bib](bibliography/biblio.bib).

Here is a simple command in order to save Matlab figures as vector graphics:

```matlab
saveas(gcf, strcat(pwd, '/img/figure'), 'epsc') % Save figure
```

The eps can then be converted to a pdf like so:

```bash
epstopdf figure.eps
```

Matlab script can also easily be added to the book:

```latex
\lstinputlisting{path/script.m}
```

You can write to me at any time either here, Facebook or via email: <lauszus@gmail.com>.
