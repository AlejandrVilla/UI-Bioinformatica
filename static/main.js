function updateDisplay() {
    const select = document.getElementById("numElements");
    const numElements = parseInt(select.value);
    const elements = document.querySelectorAll(".list-item");
    
    elements.forEach((element, index) => {
        if (index < numElements) {
            element.style.display = "block";
        } else {
            element.style.display = "none";
        }
    });
}

document.addEventListener("DOMContentLoaded", () => {
    updateDisplay();
});