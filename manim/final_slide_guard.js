(() => {
  const forwardKeys = new Set(["ArrowRight", "ArrowDown", "PageDown"]);
  const finalIndex = () => document.querySelectorAll(".slides > section").length - 1;
  const atFinalSlide = () => Reveal.getIndices().h >= finalIndex();

  function clampToFinalSlide() {
    if (Reveal.getIndices().h > finalIndex()) {
      Reveal.slide(finalIndex());
    }
  }

  Reveal.on("slidechanged", clampToFinalSlide);
  document.addEventListener("keydown", (event) => {
    if (forwardKeys.has(event.key) && atFinalSlide()) {
      event.preventDefault();
      event.stopImmediatePropagation();
    }
  }, true);
})();
